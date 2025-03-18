from docutils import nodes
from docutils.parsers.rst import directives
from sphinx.util.docutils import SphinxDirective


class TutorialCardNode(nodes.General, nodes.Element):
    """Node for a single tutorial card."""

    pass


class TutorialListNode(nodes.General, nodes.Element):
    """Node for a list of tutorial cards."""

    pass


class TutorialCardDirective(SphinxDirective):
    """
    Directive to create a tutorial card.

    Options
    -------
    title : str
        The title of the tutorial card.
    description : str
        The description of the tutorial card.
    tags : str
        The tags associated with the tutorial card.
    notebook : str
        The link to the notebook for the tutorial card.
    """

    option_spec = {
        "title": directives.unchanged,
        "description": directives.unchanged,
        "tags": directives.unchanged,
        "notebook": directives.unchanged,
    }

    def run(self):
        """
        Create a new TutorialCardNode and append it to the

        environment's all_tutorial_cards list.
        """
        env = self.state.document.settings.env
        if not hasattr(env, "all_tutorial_cards"):
            env.all_tutorial_cards = []

        title = self.options.get("title", "No Title")
        description = self.options.get("description", "")
        tags = self.options.get("tags", "")
        notebook = self.options.get("notebook", "")

        card_node = TutorialCardNode()
        card_node["title"] = title
        card_node["description"] = description
        card_node["tags"] = tags
        card_node["notebook"] = notebook

        env.all_tutorial_cards.append(card_node)
        return [card_node]


class TutorialListDirective(SphinxDirective):
    """Directive to create a list of tutorial cards."""

    def run(self):
        """Create a new TutorialListNode."""
        return [TutorialListNode("")]


def process_tutorial_cards(app, doctree, fromdocname):
    """
    Process all TutorialCardNodes and TutorialListNodes in the doctree and

    replace them with the rendered HTML.

    Parameters
    ----------
    app : Sphinx application object
        The Sphinx application object.
    doctree : docutils.nodes.document
        The doctree object.
    fromdocname : str
        The name of the document.
    """
    env = app.builder.env
    if not hasattr(env, "all_tutorial_cards"):
        return

    for node in doctree.traverse(TutorialListNode):
        container = nodes.container()

        # Render Filter Buttons
        filter_buttons = """
        <div id="filter-buttons">
            <button class="filter-btn" data-filter="all">All</button>
        """
        tags_set = set()
        for card_node in env.all_tutorial_cards:
            for tag in card_node["tags"].split(", "):
                tags_set.add(tag)
        for tag in sorted(tags_set):
            filter_buttons += f'<button class="filter-btn" data-filter="{tag}">{tag}</button>'
        filter_buttons += "</div>"

        container.append(nodes.raw("", filter_buttons, format="html"))

        # Render Tutorial Cards
        tutorial_list_html = '<div id="tutorial-cards">'
        for card_node in env.all_tutorial_cards:
            tutorial_list_html += f'''
            <div class="tutorial-card" data-tags="{card_node["tags"]}">
                <h3><a href="{card_node["notebook"]}">{card_node["title"]}</a></h3>
                <p>{card_node["description"]}</p>
                <p><strong>Tags:</strong> {card_node["tags"]}</p>
            </div>
            '''
        tutorial_list_html += "</div>"

        container.append(nodes.raw("", tutorial_list_html, format="html"))

        node.replace_self(container)


def setup(app):
    """
    App setup hook.

    Parameters
    ----------
    app : Sphinx application object
        The Sphinx application object.
    """
    app.add_node(TutorialCardNode)
    app.add_node(TutorialListNode)
    app.add_directive("tutorialcard", TutorialCardDirective)
    app.add_directive("tutoriallist", TutorialListDirective)
    app.connect("doctree-resolved", process_tutorial_cards)
