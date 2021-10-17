from input import get_tree


if __name__ == "__main__":
   tree=  get_tree(input_string="((A,B),(C,(D,E)));")
   print(tree.as_ascii_plot())
