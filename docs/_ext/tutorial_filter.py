from __future__ import annotations

from docutils import nodes

from sphinx.application import Sphinx
from sphinx.util.docutils import SphinxDirective, SphinxRole
from sphinx.util.typing import ExtensionMetadata

class TutorialListDirective(SphinxDirective):
    """A directive to create a filterable tutorial list of tutorial cards."""

    required_arguments = 1

    def run(self) -> list[nodes.Node]:
        paragraph_node = nodes.paragraph(text=f'hello {self.arguments[0]}!')
        return [paragraph_node]
    
class TutorialCardDirective(SphinxDirective):
    """A directive to create a tutorial card."""
    
    required_arguments = 1

    def run(self) -> list[nodes.Node]:
        paragraph_node = nodes.paragraph(text=f'hello {self.arguments[0]}!')
        return [paragraph_node]


def setup(app: Sphinx) -> ExtensionMetadata:
    app.add_directive("tutoriallist", TutorialListDirective)
    app.add_directive("tutorialcard", TutorialCardDirective)

    return {
        'version': '0.1',
        'parallel_read_safe': True,
        'parallel_write_safe': True,
    }