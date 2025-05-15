from docutils import nodes
from docutils.parsers.rst import Directive, directives

from sphinx.util.docutils import SphinxDirective
from docutils.statemachine import StringList

import os
import json

# Here we define the order of the model groups
# as they should appear in the list
GROUP_ORDER = [
    "Quick start",
    "scRNA-seq",
    "ATAC-seq",
    "scBS-seq",
    "Multimodal",
    "Spatial transcriptomics",
    "Model hub",
    "Common Modelling Use Cases",
    "Development",
    "R Tutorials",
]


class cardnode(nodes.container, nodes.Element):
    """A custom node representing a card in the documentation."""

    pass


class listnode(nodes.General, nodes.Element):
    """A custom node representing a list of cards in the documentation."""

    pass


def visit_card_node(self, node):
    """
    Visit a card node during the rendering process.

    Args:
        self: The visitor instance.
        node: The card node being visited.
    """
    self.visit_paragraph(node)


def depart_card_node(self, node):
    """
    Depart a card node during the rendering process.

    Args:
        self: The visitor instance.
        node: The card node being departed.
    """
    self.depart_paragraph(node)


class ListDirective(Directive):
    """A custom directive to create a list of cards in the documentation."""

    def run(self):
        """
        Process the directive and generate the list node with HTML content.

        Returns
        -------
            list: A list of nodes representing the start, placeholder, and end of the list.
        """
        list_start_html = StringList(LIST_START_HTML.split("\n"))
        start_node = nodes.paragraph()
        self.state.nested_parse(list_start_html, self.content_offset, start_node)

        list_end_html = StringList(LIST_END_HTML.split("\n"))
        end_node = nodes.paragraph()
        self.state.nested_parse(list_end_html, self.content_offset, end_node)

        # Here we insert an empty list node as a placeholder
        return start_node.children + [listnode("")] + end_node.children


class CardDirective(SphinxDirective):
    """
    A custom directive to create a card in the documentation.

    Attributes
    ----------
        has_content (bool): Indicates that the directive has content.
        option_spec (dict): Specifies the options for the directive.
    """

    has_content = True

    option_spec = {
        "path": directives.unchanged,  # link to the tutorial
        "tags": directives.unchanged,  # tags for the tutorial
    }

    def run(self):
        """
        Process the directive and generate the card node with HTML content.

        Returns
        -------
            list: A list containing the card node.
        """
        path = self.options.get("path", [])
        tags = self.options.get("tags", [])

        # Get the tutorial's title
        title = self.get_notebook_title(path)

        # Get the model group from the model group's index file
        group = self.get_index_header()

        card_node = cardnode()

        # Insert HTML content into the card node
        card_html = CARD_HTML.format(
            tags=tags,
            link=(f"{path}.html"),
            header=title,
            card_description=self.content[0],
            group=group,
        )

        card_list = StringList(card_html.split("\n"))
        self.state.nested_parse(card_list, self.content_offset, card_node)

        env = self.state.document.settings.env

        if not hasattr(env, "all_cards"):
            env.all_cards = []

        env.all_cards.append(card_node)

        return [card_node]

    def get_notebook_title(self, path):
        """Retrieve the title of the notebook directly from the `.ipynb` file."""
        # Determine the actual file path (before processing)
        docs_root = self.env.srcdir  # Root directory of the source docs
        tutorials_dir = os.path.join(docs_root, "tutorials/")
        notebook_path = os.path.join(tutorials_dir, path)

        # Ensure the path has the `.ipynb` extension
        if not notebook_path.endswith(".ipynb"):
            notebook_path += ".ipynb"

        if not os.path.exists(notebook_path):
            print(f"Notebook file {notebook_path} not found.")
            return f"Notebook not found, path tried: {notebook_path}"

        # Read the JSON content of the notebook
        try:
            with open(notebook_path, encoding="utf-8") as f:
                notebook_data = json.load(f)

            # Extract the first-level heading from the notebook cells
            for cell in notebook_data.get("cells", []):
                if cell.get("cell_type") == "markdown":
                    for line in cell.get("source", []):
                        if line.startswith("# "):  # First-level heading
                            return line.lstrip("# ").strip()
        except Exception as e:
            print(f"Error reading notebook file {notebook_path}: {e}")
            return f"Error reading notebook: {e}"

        return "No title found"

    def get_index_header(self):
        """
        Retrieve the title of the model group index from its top-level header.

        Returns
        -------
            str: The title of the model group index, or "No title found" if unavailable.
        """
        # Get the title from the model group index's top level header (such as scATAC-seq)
        doctree = self.state.document

        # get the first title node
        for node in doctree.traverse(nodes.title):
            return node.astext()

        return "No title found"


def process_card_nodes(app, doctree, fromdocname):
    """
    Process card nodes and replace them with their rendered content in the list.

    Args:
        app: The Sphinx application instance.
        doctree: The document tree being processed.
        fromdocname: The name of the document being processed.
    """
    env = app.builder.env

    if not hasattr(env, "all_cards"):
        env.all_cards = []

    # Sort the cards by group name based on the GROUP_ORDER
    def group_sort_key(card):
        # Find group name (header of the page)
        group = "No group found"
        for node in doctree.traverse(nodes.title):
            group = node.astext()

        group_index = GROUP_ORDER.index(group) if group in GROUP_ORDER else len(GROUP_ORDER)
        return (group_index, card.source, card.line)

    # Sort the cards by group name
    env.all_cards.sort(key=group_sort_key)

    # Don't want card to render where the directive is, but in the list
    for card in doctree.traverse(cardnode):
        card.parent.remove(card)

    # Loop through all list directive nodes. Only 1 for now (for tutorials)
    for node in doctree.traverse(listnode):
        content = []

        for card in env.all_cards:
            content.append(card)

        node.replace_self(content)


def setup(app):
    """
    Set up the Sphinx extension by adding nodes, directives, and event handlers.

    Args:
        app: The Sphinx application instance.

    Returns
    -------
        dict: Metadata about the extension.
    """
    app.add_node(listnode)
    app.add_node(cardnode, html=(visit_card_node, depart_card_node))
    app.add_directive("customlist", ListDirective)
    app.add_directive("customcard", CardDirective)
    app.connect("doctree-resolved", process_card_nodes)

    return {
        "parallel_read_safe": True,
        "parallel_write_safe": True,
    }


# html templates for list and cards

CARD_HTML = """
<div class="card-container" data-tags="{tags}" data-group="{group}">
    <div class="card">
        <a href="{link}">
            <div class="card-body">
                <div class="card-title-container">
                    <h4>{header}</h4>
                </div>
                <p class="card-description">{card_description}</p>
                <p class="model-group-name">{group}</p>
            </div>
        </a>
    </div>
</div>
"""


LIST_START_HTML = """
<div id="cards-container">

<nav class="navbar tutorials-nav">
    <div class="tags-container">
        <div id="dropdown-filter-tags">
            <div class="filter-menu all-tag-selected">
                <div class="filter filter-btn" data-tag="all">All</div>
            </div>
        </div>
    </div>
</nav>

<hr class="tutorials-hr">

<nav class="navbar">
    <div class="tabs-container">
        <div class="tab-menu">
            <div class="tab tab-selected" data-group="all">All</div>
        </div>
    </div>
</nav>

<div class="row">

<div id="cards">
<div class="list">
"""

LIST_END_HTML = """
</div>

<div class="pagination d-flex justify-content-center"></div>

</div>

</div>
<br>
<br>
"""
