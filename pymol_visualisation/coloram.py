from pymol import cmd
import colorsys


def interpolate_color(value, min_val=0.0, max_val=1.0):
    """
    Interpolates between blue (0), grey (0.5), and red (1.0).
    """
    # Clamp value
    v = max(min(value, max_val), min_val)

    # Colors in RGB
    blue = [33 / 255.0, 102 / 255.0, 172 / 255.0]
    grey = [230 / 255.0, 230 / 255.0, 230 / 255.0]
    red = [178 / 255.0, 24 / 255.0, 43 / 255.0]

    if v <= 0.5:
        # Interpolate between blue and grey
        t = v / 0.5
        color = [blue[i] * (1 - t) + grey[i] * t for i in range(3)]
    else:
        # Interpolate between grey and red
        t = (v - 0.5) / 0.5
        color = [grey[i] * (1 - t) + red[i] * t for i in range(3)]

    return color


def coloram(selection="all"):
    """
    Colors structures using a gradient based on B-factor.
    """
    model = cmd.get_model(selection)
    for i, atom in enumerate(model.atom):
        b = atom.b
        color = interpolate_color(b)
        cname = f"amcolor_{i}"
        cmd.set_color(cname, color)
        cmd.color(cname, f"id {atom.index}")


cmd.extend("coloram", coloram)
cmd.auto_arg[0]["coloram"] = [cmd.object_sc, "object", ""]
