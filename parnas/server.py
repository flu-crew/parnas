# -*- coding: utf-8 -*-
import io
import os
import re
import sys
import random as rnd
import tempfile
import traceback
from math import floor

try:
    from flask import Flask, request, jsonify, send_from_directory
except ImportError:
    raise ImportError(
        "Flask is required for PARNAS server mode.\n"
        "Install it with: pip install flask"
    )

from dendropy import Tree
from scipy.stats import percentileofscore

from parnas.cli import color_by_clusters
from parnas.medoids import (
    build_distance_functions,
    binarize_tree,
    get_costs,
    find_n_medoids,
    find_n_medoids_with_diversity,
    find_coverage,
)
from parnas.medoids.medoid_utils import get_centers_score, compute_percent_coverage
from parnas.options import reweigh_tree_ancestral

# When frozen by PyInstaller (--onefile), __file__ is the executable itself
# and all package data is extracted to sys._MEIPASS.
if getattr(sys, "frozen", False):
    WEB_DIR = os.path.join(sys._MEIPASS, "parnas", "web")  # type: ignore[attr-defined]
else:
    WEB_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "web")
app = Flask(__name__)


# ── Helpers ───────────────────────────────────────────────────────────────────

def _match_taxa(tree: Tree, regex: str):
    """Return list of taxon labels matching regex; [] if empty; raises ValueError on bad regex."""
    if not regex:
        return []
    try:
        regex_compiled = re.compile(regex)
        return [t.label for t in tree.taxon_namespace if regex_compiled.fullmatch(t.label)]
    except re.error as exc:
        raise ValueError(f"Invalid regex '{regex}': {exc}") from exc


def _parse_weights(path: str) -> dict:
    """Parse a weights CSV; raises ValueError on malformed input."""
    weights = {}
    with open(path) as f:
        headers = [h.strip() for h in f.readline().split(",")]
        if headers[:2] != ["taxon", "weight"]:
            raise ValueError("Weights CSV must have columns 'taxon' and 'weight'")
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = [p.strip() for p in line.split(",")]
            if len(parts) < 2:
                raise ValueError(f"Malformed weights row: {line!r}")
            taxon, weight_str = parts[0], parts[1]
            try:
                w = float(weight_str)
            except ValueError:
                raise ValueError(f"Non-numeric weight '{weight_str}' for taxon '{taxon}'")
            if w <= 0 or w >= 1000:
                raise ValueError(f"Weight {w} for '{taxon}' out of range [0, 1000]")
            weights[taxon] = w
    return weights


def _extract_colors(tree: Tree) -> dict:
    return {
        t.label: t.annotations.get_value("!color")
        for t in tree.taxon_namespace
        if t.annotations.get_value("!color")
    }


def _newick(tree: Tree) -> str:
    buf = io.StringIO()
    tree.write(file=buf, schema="newick")
    return buf.getvalue().strip()


# ── Routes ────────────────────────────────────────────────────────────────────

@app.route("/")
def index():
    return send_from_directory(WEB_DIR, "index.html")


@app.route("/api/run", methods=["POST"])
def run_parnas():
    tree_file = request.files.get("tree")
    if not tree_file:
        return jsonify({"error": "No tree file provided"}), 400

    # ── Boolean flags ─────────────────────────────────────────
    cover     = request.form.get("cover")    == "true"
    evaluate  = request.form.get("evaluate") == "true"
    is_binary = request.form.get("binary")   == "true"

    # ── Regex strings ─────────────────────────────────────────
    prior_regex     = request.form.get("prior",           "").strip()
    exc_rep_regex   = request.form.get("exclude_rep",     "").strip()
    exc_obj_regex   = request.form.get("exclude_obj",     "").strip()
    exc_full_regex  = request.form.get("exclude_fully",   "").strip()
    constrain_regex = request.form.get("constrain_fully", "").strip()

    # ── Numeric params ────────────────────────────────────────
    radius_str    = request.form.get("radius",    "").strip()
    threshold_str = request.form.get("threshold", "").strip()

    n = -1
    if not cover and not evaluate:
        try:
            n = int(request.form.get("n", ""))
        except (ValueError, TypeError):
            return jsonify({"error": "n must be a positive integer (or enable --cover)"}), 400

    radius = None
    if radius_str:
        try:
            radius = float(radius_str)
            if radius <= 0:
                return jsonify({"error": "Radius must be a positive number"}), 400
        except ValueError:
            return jsonify({"error": "Radius must be a valid number"}), 400

    threshold = None
    if threshold_str:
        try:
            threshold = float(threshold_str)
            if threshold <= 0 or threshold >= 100:
                return jsonify({"error": "Threshold must be between 0 and 100 (exclusive)"}), 400
        except ValueError:
            return jsonify({"error": "Threshold must be a valid number"}), 400

    if cover and radius is None and threshold is None:
        return jsonify({"error": "--cover requires --radius or --threshold"}), 400

    # ── Temp file management ──────────────────────────────────
    temp_files = []

    def save_upload(file_obj, ext=".tmp"):
        fd, path = tempfile.mkstemp(suffix=ext)
        os.close(fd)
        file_obj.save(path)
        temp_files.append(path)
        return path

    tree_path = save_upload(tree_file, os.path.splitext(tree_file.filename)[1] or ".tre")

    try:
        # ── Load tree ─────────────────────────────────────────
        tree = None
        last_err = None
        for schema in ("newick", "nexus"):
            try:
                tree = Tree.get(path=tree_path, schema=schema, preserve_underscores=True)
                break
            except Exception as e:
                last_err = e
        if tree is None:
            return jsonify({"error": f"Cannot parse tree file: {last_err}"}), 400

        n_taxa = len(tree.taxon_namespace)
        if not cover and not evaluate and (n < 1 or n >= n_taxa):
            return jsonify({
                "error": f"n must be between 1 and {n_taxa - 1} (tree has {n_taxa} taxa)"
            }), 400

        # ── Resolve regex fields ──────────────────────────────
        try:
            prior_centers  = _match_taxa(tree, prior_regex) or None
            excluded_taxa  = _match_taxa(tree, exc_rep_regex)
            obj_excluded   = _match_taxa(tree, exc_obj_regex)
            fully_excluded = _match_taxa(tree, exc_full_regex)
            if constrain_regex:
                constrained = _match_taxa(tree, constrain_regex)
                all_labels  = {t.label for t in tree.taxon_namespace}
                fully_excluded += list(all_labels - set(constrained))
        except ValueError as exc:
            return jsonify({"error": str(exc)}), 400

        if evaluate and not prior_centers:
            return jsonify({
                "error": "--evaluate requires at least one prior representative (check --prior regex)"
            }), 400

        # ── Weights ───────────────────────────────────────────
        taxa_weights = None
        weights_file = request.files.get("weights")
        if weights_file and weights_file.filename:
            w_path = save_upload(weights_file, ".csv")
            try:
                taxa_weights = _parse_weights(w_path)
            except ValueError as exc:
                return jsonify({"error": f"Weights file error: {exc}"}), 400

        # ── Threshold / tree re-weighting ─────────────────────
        query_tree = tree
        if threshold is not None:
            aln_type = request.form.get("aln_type", "nt")
            aln_file = request.files.get("alignment")
            if not aln_file or not aln_file.filename:
                return jsonify({"error": "--threshold requires an alignment file"}), 400
            is_aa  = aln_type == "aa"
            a_path = save_upload(aln_file, os.path.splitext(aln_file.filename)[1] or ".fasta")
            try:
                from Bio import AlignIO
                alignment  = list(AlignIO.read(a_path, "fasta"))
                radius     = floor((1 - threshold / 100) * len(alignment[0]))
                query_tree = reweigh_tree_ancestral(tree_path, a_path, is_aa)
            except Exception as exc:
                return jsonify({"error": f"Threshold/TreeTime error: {exc}"}), 400

        # ── Common setup ──────────────────────────────────────
        binarize_tree(query_tree, edge_length=0)
        cost_map           = get_costs(query_tree, excluded_taxa, fully_excluded)
        fully_excluded_all = list(set(fully_excluded) | set(obj_excluded))

        # ── Evaluate mode ─────────────────────────────────────
        if evaluate:
            if cover:
                pct = compute_percent_coverage(
                    query_tree, prior_centers, radius, fully_excluded_all
                )
                color_by_clusters(
                    query_tree, prior_centers, prior_centers=[],
                    fully_excluded=fully_excluded, radius=radius,
                )
                return jsonify({
                    "mode":          "evaluate_cover",
                    "tree":          _newick(query_tree),
                    "prior_centers": prior_centers,
                    "colors":        _extract_colors(query_tree),
                    "coverage_pct":  round(pct * 100, 2),
                    "n_taxa":        n_taxa,
                })

            dist_fns = build_distance_functions(
                query_tree, is_binary=is_binary,
                fully_excluded=fully_excluded_all, radius=radius,
                taxa_weights=taxa_weights,
            )
            prior_score = get_centers_score(query_tree, prior_centers, dist_fns)
            n_eval      = len(prior_centers)
            reps, _, div_scores, obj1 = find_n_medoids_with_diversity(
                query_tree, n_eval, dist_fns, cost_map, max_dist=radius
            )

            prior_diversity = (obj1 - prior_score) / obj1 * 100 if obj1 > 0 and n_eval > 1 else None
            best_diversity  = div_scores[-1] if div_scores and n_eval > 1 else None

            taxa_labels = [
                lf.taxon.label for lf in query_tree.leaf_nodes()
                if lf.taxon.label not in fully_excluded_all
            ]
            rnd_scores = []
            for _ in range(1000):
                rnd.shuffle(taxa_labels)
                rnd_scores.append(get_centers_score(query_tree, taxa_labels[:n_eval], dist_fns))
            percentile = 100 - percentileofscore(sorted(rnd_scores), prior_score, kind="strict")

            color_by_clusters(
                query_tree, reps, prior_centers=prior_centers,
                fully_excluded=fully_excluded, radius=radius,
            )
            return jsonify({
                "mode":                  "evaluate",
                "tree":                  _newick(query_tree),
                "prior_centers":         prior_centers,
                "best_representatives":  reps,
                "colors":                _extract_colors(query_tree),
                "prior_diversity":       round(prior_diversity, 2) if prior_diversity is not None else None,
                "best_diversity":        round(best_diversity,  2) if best_diversity  is not None else None,
                "percentile":            round(percentile, 1),
                "n_taxa":                n_taxa,
            })

        # ── Sample mode ───────────────────────────────────────
        dist_fns = build_distance_functions(
            query_tree, prior_centers=prior_centers, is_binary=is_binary,
            fully_excluded=fully_excluded_all, radius=radius,
            taxa_weights=taxa_weights,
        )

        if cover:
            coverage = find_coverage(
                query_tree, radius, cost_map, prior_centers, fully_excluded, obj_excluded
            )
            if coverage is None:
                reps, prev_value = [], -1
                for k in range(1, n_taxa + 1):
                    reps_k, value = find_n_medoids(query_tree, k, dist_fns, cost_map, max_dist=radius)
                    if value == 0:
                        reps = reps_k
                        break
                    if value == prev_value:
                        break
                    reps, prev_value = reps_k, value
                representatives = reps
            else:
                representatives = coverage
            diversity_scores = None
        else:
            representatives, _, diversity_scores, _ = find_n_medoids_with_diversity(
                query_tree, n, dist_fns, cost_map, max_dist=radius
            )

        color_by_clusters(
            query_tree, representatives,
            prior_centers=prior_centers,
            fully_excluded=fully_excluded,
            radius=radius,
        )

        clusters = {
            lf.taxon.label: int(lf.annotations.get_value("center"))
            for lf in query_tree.leaf_nodes()
            if lf.annotations.get_value("center") is not None
            and int(lf.annotations.get_value("center")) >= 0
        }

        return jsonify({
            "mode":            "sample",
            "representatives": representatives,
            "prior_centers":   prior_centers or [],
            "colors":          _extract_colors(query_tree),
            "clusters":        clusters,
            "tree":            _newick(query_tree),
            "diversity":       float(diversity_scores[-1]) if diversity_scores else None,
            "n_taxa":          n_taxa,
        })

    except Exception as exc:
        return jsonify({"error": str(exc), "traceback": traceback.format_exc()}), 500
    finally:
        for p in temp_files:
            try:
                os.unlink(p)
            except OSError:
                pass


def run_server(host: str = "localhost", port: int = 8080) -> None:
    print(f"\n  PARNAS web server  →  http://{host}:{port}\n")
    app.run(host=host, port=port, debug=False)
