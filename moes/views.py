import sys
from django.http import HttpResponse
from django.template import loader

sys.path.append("src/")

from src.specdata import SpecDB


def index(request):
    template = loader.get_template("moes/index.html")
    context = {}

    if request.method == "POST":
        trot_input = int(request.POST.get("Trot_input"))
        spec = get_spectrum(trot_input)
        context["result"] = spec
    return HttpResponse(template.render(context, request))


def get_spectrum(Trot):
    OHAX = SpecDB("OHAX.db")
    Tvib = 1000
    wmin = 310
    wmax = 320
    spec = OHAX.get_spectrum(Trot=Trot, Tvib=Tvib, wmin=wmin, wmax=wmax)
    spec.refine_mesh(points_per_nm=50)
    spec.convolve_with_slit_function(
        gauss=0.0396, lorentz=0.0138, instrumental_step=0.0318
    )
    with open("moes.log", "w") as f:
        f.write(f"X values: {spec.x}\n")

    # return {"x": spec.x.tolist(), "y": spec.y.tolist()}

    return [{"x": x, "y": y} for x, y in zip(spec.x.tolist(), spec.y.tolist())]
