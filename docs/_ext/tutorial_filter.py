from docutils import nodes
from docutils.statemachine import StringList

from sphinx.application import Sphinx
from sphinx.util.docutils import SphinxDirective
from sphinx.util.typing import ExtensionMetadata


class tutorialcardnode(nodes.General, nodes.TextElement):
    """Node for a single tutorial card."""

    pass


class tutoriallistnode(nodes.General, nodes.Element):
    """Node for the list of tutorial cards."""

    pass


def visit_tutorialcard_node(self, node):
    """Visit a tutorial card node."""
    self.visit_paragraph(node)


def depart_tutorialcard_node(self, node):
    """Depart from a tutorial card node."""
    self.depart_paragraph(node)


def purge_tutorial_cards(app, env, docname):
    """Purge tutorial cards from the environment when a document is removed."""
    if not hasattr(env, "tutorial_card_node_list"):
        return

    env.tutorial_card_node_list = [
        node for node in env.tutorial_card_node_list if node["docname"] != docname
    ]


def merge_tutorial_cards(app, env, docnames, other):
    """Merge tutorial cards from different environments."""
    if not hasattr(env, "tutorial_card_node_lists"):
        env.tutorial_card_node_lists = []
    if hasattr(other, "tutorial_card_node_list"):
        env.tutorial_card_node_list.extend(other.tutorial_card_node_list)


class TutorialListDirective(SphinxDirective):
    """Directive to insert a list of tutorial cards."""

    def run(self):
        """Run the directive to insert header and footer markup for tutorial card list."""
        header_list = StringList(TUTORIAL_LIST_START.split("\n"))
        header_node = nodes.paragraph()
        self.state.nested_parse(header_list, self.content_offset, header_node)

        footer_list = StringList(TUTORIAL_LIST_END.split("\n"))
        footer_node = nodes.paragraph()
        self.state.nested_parse(footer_list, self.content_offset, footer_node)

        return header_node.children + [tutoriallistnode("")] + footer_node.children


# TODO: HTML templates below borrowed from ____

# HTML template for tutorial list header
TUTORIAL_LIST_START = """
.. raw:: html

   <div id="tutorial-cards-container">

   <nav class="navbar navbar-expand-lg navbar-light tutorials-nav col-12">
     <div class="tutorial-tags-container">
         <div id="dropdown-filter-tags">
             <div class="tutorial-filter-menu">
                 <div class="tutorial-filter filter-btn all-tag-selected" data-tag="all">All</div>
             </div>
         </div>
     </div>
   </nav>

   <hr class="tutorials-hr">

   <div class="row">

   <div id="tutorial-cards">
   <div class="list">

.. Tutorial cards below this line

"""

# HTML template for tutorial list footer
TUTORIAL_LIST_END = """

.. End of tutorial card section

.. raw:: html

   </div>

   <div class="pagination d-flex justify-content-center"></div>

   </div>

   </div>
   <br>
   <br>
"""


class TutorialCardDirective(SphinxDirective):
    """Directive to create a single tutorial card."""

    option_spec = {
        "link": str,
        "description": str,
        "tags": str,
    }

    def run(self) -> list[nodes.Node]:
        """Run the directive to create a tutorial card."""
        link = self.options["link"]
        description = self.options["description"]

        if "tags" in self.options:  # TODO: how should these be processed?
            tags = self.options["tags"]

        doctree = self.env.get_doctree(link)

        # Try to get the title from the environment's titles dictionary
        title_node = self.env.titles.get(link)
        if title_node:
            title = title_node.astext()  # Get the text of the title node
        else:
            title = "Untitled"  # Fallback if the title isn't found

        # Get the model group (scRNA-seq, scATAC-seq, etc.)
        doctree = self.env.get_doctree(self.env.docname)
        for section in doctree.traverse(nodes.section):
            model_group = section[0].astext()
            break

        # Create all nodes for tutorial card
        link = f"{link}.html"

        tutorial_card_rst = TUTORIAL_CARD_TEMPLATE.format(
            link=link,
            header=title,
            card_description=description,
            tags=tags,
            model_group=model_group,
        )
        tutorial_card_list = StringList(tutorial_card_rst.split("\n"))
        node = tutorialcardnode()
        self.state.nested_parse(tutorial_card_list, self.content_offset, node)

        if not hasattr(self.env, "tutorial_card_node_list"):
            self.env.tutorial_card_node_list = []
        self.env.tutorial_card_node_list.append({"docname": self.env.docname, "node": node})
        return [node]


# HTML template for each tutorial card
TUTORIAL_CARD_TEMPLATE = """
.. raw:: html

    <div class="tutorial-card-container">
        <div class="tutorial-card">
            <a href="{link}" style="text-decoration: none;">
                <div class="tutorial-card-body">
                    <h4 class="tutorial-card-title">{header} {beta}</h4>
                    <p class="tutorial-card-description">{card_description}</p>
                </div>
            </a>
        </div>
    </div>
"""


def process_tutorial_card_nodes(app, doctree, fromdocname):
    """Process tutorial card nodes in the document tree."""
    # Remove all tutorial card nodes in the document. We don't want these to display
    for card in doctree.traverse(tutorialcardnode):
        card.parent.remove(card)

    # Add the actual tutorial cards to the tutorial list node
    for node in doctree.traverse(tutoriallistnode):
        content = []

        for card in app.builder.env.tutorial_card_node_list:
            content.append(card["node"])

        node.replace_self(content)


def setup(app: Sphinx) -> ExtensionMetadata:
    """Set up the Sphinx extension."""
    app.add_node(tutoriallistnode)
    app.add_node(
        tutorialcardnode,
        html=(visit_tutorialcard_node, visit_tutorialcard_node),
        latex=(visit_tutorialcard_node, visit_tutorialcard_node),
        text=(visit_tutorialcard_node, visit_tutorialcard_node),
    )

    app.add_directive("tutorialcard", TutorialCardDirective)
    app.add_directive("tutoriallist", TutorialListDirective)

    app.connect("env-purge-doc", purge_tutorial_cards)
    app.connect("env-merge-info", merge_tutorial_cards)
    app.connect("doctree-resolved", process_tutorial_card_nodes)

    return {
        "parallel_read_safe": True,
        "parallel_write_safe": True,
    }
