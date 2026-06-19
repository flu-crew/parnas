/**
 * Newick/NEXUS parser → immutable node tree.
 * Node shape: { id, name, length, support, children, annotations }
 * annotations: object parsed from [&key=value,...] NEXUS comments
 */

let _idCounter = 0;

function makeNode(name = "", length = null, support = null) {
  return { id: _idCounter++, name, length, support, children: [], annotations: {} };
}

/** Parse [&key=value,...] annotation comment string (without brackets). */
function parseAnnotations(raw) {
  if (!raw || !raw.startsWith("&")) return {};
  const out = {};
  const inner = raw.slice(1);
  // Simple key=value pairs, values may be quoted or unquoted
  const re = /(\w+)=("([^"]*?)"|([^,}\]]+))/g;
  let m;
  while ((m = re.exec(inner)) !== null) {
    out[m[1]] = m[3] !== undefined ? m[3] : m[4];
  }
  return out;
}

/**
 * Tokenise a Newick string into tokens:
 * '(' ')' ',' ';' label branch_length annotation
 */
function tokenise(s) {
  const tokens = [];
  let i = 0;
  while (i < s.length) {
    const ch = s[i];
    if (ch === "(" || ch === ")" || ch === "," || ch === ";") {
      tokens.push({ type: ch, val: ch });
      i++;
    } else if (ch === ":") {
      // branch length
      i++;
      let num = "";
      while (i < s.length && /[0-9eE+\-.]/.test(s[i])) num += s[i++];
      tokens.push({ type: "length", val: parseFloat(num) });
    } else if (ch === "[") {
      // annotation comment
      i++;
      let inner = "";
      while (i < s.length && s[i] !== "]") inner += s[i++];
      i++; // skip ]
      tokens.push({ type: "annotation", val: inner });
    } else if (ch === "'") {
      // single-quoted label
      i++;
      let label = "";
      while (i < s.length) {
        if (s[i] === "'" && s[i + 1] === "'") { label += "'"; i += 2; }
        else if (s[i] === "'") { i++; break; }
        else label += s[i++];
      }
      tokens.push({ type: "label", val: label });
    } else if (/\s/.test(ch)) {
      i++;
    } else {
      // unquoted label / internal node name
      let label = "";
      while (i < s.length && !/[(),;:\[\]'"]/.test(s[i]) && !/\s/.test(s[i])) {
        label += s[i++];
      }
      if (label) tokens.push({ type: "label", val: label });
    }
  }
  return tokens;
}

/**
 * Parse a Newick string, return root node.
 * Raises on malformed input.
 */
export function parseNewick(newick) {
  _idCounter = 0;
  // Strip NEXUS wrapper if present
  const treeMatch = newick.match(/tree\s+\S+\s*=\s*(\[.*?\]\s*)?([^;]+;)/i);
  const src = (treeMatch ? treeMatch[2] : newick).trim();

  const tokens = tokenise(src);
  let pos = 0;

  function peek() { return tokens[pos]; }
  function consume(type) {
    const t = tokens[pos++];
    if (type && t?.type !== type) throw new Error(`Expected ${type}, got ${t?.type} at token ${pos}`);
    return t;
  }

  function parseNode() {
    const node = makeNode();
    if (peek()?.type === "(") {
      consume("(");
      node.children.push(parseNode());
      while (peek()?.type === ",") {
        consume(",");
        node.children.push(parseNode());
      }
      consume(")");
    }
    // label / support after closing paren or leaf
    if (peek()?.type === "annotation") {
      node.annotations = parseAnnotations(consume("annotation").val);
    }
    if (peek()?.type === "label") {
      const lbl = consume("label").val;
      if (node.children.length > 0) {
        // internal node: label may be support value
        const num = parseFloat(lbl);
        if (!isNaN(num) && String(num) === lbl) node.support = num;
        else node.name = lbl;
      } else {
        node.name = lbl;
      }
    }
    if (peek()?.type === "annotation") {
      Object.assign(node.annotations, parseAnnotations(consume("annotation").val));
    }
    if (peek()?.type === "length") {
      node.length = consume("length").val;
    }
    if (peek()?.type === "annotation") {
      Object.assign(node.annotations, parseAnnotations(consume("annotation").val));
    }
    return node;
  }

  const root = parseNode();
  return root;
}

/**
 * Serialise to NEXUS format, embedding per-leaf annotation blocks.
 * extraAnnotations: Map<leafName, {key: value}>
 */
export function toNexus(root, extraAnnotations = new Map()) {
  function annotatedNewick(node) {
    let s = "";
    if (node.children.length > 0) {
      s += "(" + node.children.map(annotatedNewick).join(",") + ")";
    }
    // Merge stored annotations with any extra ones keyed by leaf name
    const ann = { ...(node.annotations || {}) };
    if (node.children.length === 0 && node.name && extraAnnotations.has(node.name)) {
      Object.assign(ann, extraAnnotations.get(node.name));
    }
    const annKeys = Object.keys(ann);
    if (annKeys.length > 0) {
      s += `[&${annKeys.map(k => `${k}=${ann[k]}`).join(",")}]`;
    }
    if (node.name) s += node.name.includes(" ") ? `'${node.name}'` : node.name;
    if (node.support !== null && node.support !== undefined) s += node.support;
    if (node.length !== null && node.length !== undefined) s += `:${node.length}`;
    return s;
  }

  const nwk = annotatedNewick(root);
  return `#NEXUS\nBegin trees;\n\ttree parnas = [&R] ${nwk};\nEnd;\n`;
}
