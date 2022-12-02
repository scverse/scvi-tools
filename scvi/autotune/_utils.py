import sys


def in_notebook() -> bool:
    """
    Check if running in a Jupyter notebook or Colab session.

    Based on: https://stackoverflow.com/questions/15411967/how-can-i-check-if-code-is-executed-in-the-ipython-notebook
    """
    try:
        from IPython import get_ipython

        shell = get_ipython().__class__.__name__
        if shell == "ZMQInteractiveShell":
            return True
        elif shell == "TerminalInteractiveShell":
            return False
        else:
            return False
    except ImportError:
        in_colab = "google.colab" in sys.modules
        return in_colab
