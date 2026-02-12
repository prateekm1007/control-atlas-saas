GLOBAL_WINDOWS = {("C","C"):(1.2,1.75),("C","N"):(1.2,1.65),("C","O"):(1.1,1.6),("C","S"):(1.6,1.95),("N","O"):(1.1,1.55),("P","O"):(1.4,1.75)}
CONTEXT_WINDOWS = {
    ("S","S"): {"allowed":["CYS"], "win":(1.85,2.25)},
    ("ZN","N"): {"allowed":["HIS","CYS"], "win":(1.8,2.35)},
    ("ZN","O"): {"allowed":["ASP","GLU","HOH"], "win":(1.8,2.35)},
    ("MG","O"): {"allowed":["ASP","GLU","HOH","ATP","ADP"], "win":(1.8,2.45)}
}
def is_contextually_bonded(a1, a2, d):
    p = tuple(sorted([a1.element.upper(), a2.element.upper()]))
    if p in CONTEXT_WINDOWS:
        r = CONTEXT_WINDOWS[p]
        if a1.res_name.upper() in r["allowed"] or a2.res_name.upper() in r["allowed"]:
            return r["win"][0] <= d <= r["win"][1]
    return GLOBAL_WINDOWS[p][0] <= d <= GLOBAL_WINDOWS[p][1] if p in GLOBAL_WINDOWS else False
