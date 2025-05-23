import io
import base64
import matplotlib.pyplot as plt

def _fig_to_base64(obj) -> str:

    if isinstance(obj, (list, tuple)) and obj:
        obj = obj[0]

    if hasattr(obj, "savefig"):
        fig = obj
    elif hasattr(obj, "figure"):
        fig = obj.figure
    elif hasattr(obj, "fig"):
        fig = obj.fig
    else:
        raise TypeError(f"Cannot extract Figure from object of type {type(obj)}")

    buf = io.BytesIO()
    fig.savefig(buf, format="png", bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    return base64.b64encode(buf.read()).decode("ascii")
