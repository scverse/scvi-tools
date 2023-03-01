from sphinx.application import Sphinx


def edit_colab_url(
    app: Sphinx,
    pagename: str,
    templatename: str,
    context: dict,
    doctree: str,
):
    """Edit the colab url to point to the correct repo.

    This assumes that the tutorials repo makes the same tag releases as the main repo.
    """
    try:
        header_buttons = context["header_buttons"]
    except KeyError:
        return
    for button in header_buttons:
        # get launch buttons
        if button["label"] == "launch-buttons":
            # only one items in the launch buttons list as we only use colab
            # remove "tutorials/notebooks" from url
            button["buttons"][0]["url"] = button["buttons"][0]["url"].replace(
                "/docs/tutorials/notebooks", ""
            )
            button["buttons"][0]["url"] = button["buttons"][0]["url"].replace(
                "scvi-tools", "scvi-tutorials"
            )


def setup(app: Sphinx):
    """Setup the extension."""
    app.connect("html-page-context", edit_colab_url)
