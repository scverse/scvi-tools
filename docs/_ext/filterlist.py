from docutils import nodes
from docutils.parsers.rst import Directive, directives

from sphinx.util.docutils import SphinxDirective
from docutils.statemachine import StringList

import os


class cardnode(nodes.Container, nodes.Element):
    """A custom node representing a card in the documentation."""


class listnode(nodes.General, nodes.Element):
    """A custom node representing a list of cards in the documentation."""


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
        # tags = self.options.get("tags", [])

        # Get the tutorial's title
        title = self.get_notebook_title(path)

        # Get the model group from the model group's index file
        group = self.get_index_header()

        card_node = cardnode()

        # Insert HTML content into the card node
        card_html = CARD_HTML.format(
            link=path, header=title, card_description=self.content[0], group=group
        )

        card_list = StringList(card_html.split("\n"))
        self.state.nested_parse(card_list, self.content_offset, card_node)

        env = self.state.document.settings.env
        if not hasattr(env, "all_cards"):
            env.all_cards = []

        env.all_cards.append(card_node)

        return [card_node]

    def get_notebook_title(self, path):
        """
        Retrieve the title of the notebook from its top-level header.

        Args:
            path (str): The path to the notebook.

        Returns
        -------
            str: The title of the notebook, or "No title found" if unavailable.
        """
        # First get parent directory of where the card is, so we can get the full path
        doc_path = self.env.docname
        parent_dir = os.path.dirname(doc_path)

        # Get the title from the notebook's top level header
        doctree = self.env.get_doctree(f"{parent_dir}/{path}")

        # get the first title node
        for node in doctree.traverse(nodes.title):
            return node.astext()

        return "No title found"

    def get_index_header(self):
        """
        Retrieve the title of the model group index from its top-level header.

        Returns
        -------
            str: The title of the model group index, or "No title found" if unavailable.
        """
        # Get the title from the model group index's top level header (such as scATAC-seq)
        doctree = self.env.get_doctree(self.env.docname)

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
<div class="card-container">
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
            <div class="filter-menu">
                <div class="filter filter-btn all-tag-selected" data-tag="all">All</div>
            </div>
        </div>
    </div>
</nav>

<hr class="tutorials-hr">

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
