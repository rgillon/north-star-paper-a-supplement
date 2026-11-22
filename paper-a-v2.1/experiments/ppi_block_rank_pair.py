"""Verify the donor systems and exact sections for blockwise rank interpolation."""
from fractions import Fraction
from math import gcd, prod

from ppi_coarse_pair import select_sections
from ppi_rank_box_pair import blend, odd_part
from ppi_reallocated_donors import packet_bound
from ppi_reduced_lower import MASKS_1, MASKS_2, SIGNED_1, SIGNED_2, PATH, pinned

UPPER_DONORS = [
    ((3, 5, 7, 11, 17, 53),
     [[113, 9799832791439900672], [233, 18432257494765043712], [163, 12297826453260599296], [97, 9223372039002259456], [191, 16717537740941426688], [127, 11529355786241212416], [103, 9223512776490647552], [229, 18442240478376099840], [239, 18446649518561724416], [-197, 17870283325566877696], [-151, 13835198795991777280], [-149, 12150711797474656256], [-139, 11574391782525403136], [-211, 18374686483949813760], [-137, 12105675801190465536], [-241, 18446687451151370368], [-223, 18419722480234004480], [-199, 17361376567555522560]],
     [2563, 4674, 4872, 5156, 6168, 10274, 17432, 20485, 33924, 66192, 81954, 131590, 147522, 7305, 7496, 75011, 163945],
     4240, 199),
    ((13, 19, 29, 31, 47),
     [[521641, 2147483648], [540271, 2281701376], [550621, 2684354560], [561961, 2684387328], [562591, 2684395520], [565111, 3221225472], [568999, 2290122752], [569251, 2818572288], [570403, 2290649216], [580357, 2692784256], [580693, 2694873216], [581323, 2694881408], [582499, 3229614080], [594091, 3758096384], [596611, 3233857664], [598987, 2827001984], [652063, 4293847240], [650479, 4243124384], [651331, 4276157120], [-651401, 4277205664], [-651461, 4294495432], [-651509, 4293983424], [-651587, 4277854400], [-651697, 4293849280], [-651767, 4294897824], [-651877, 4276676800], [-651901, 4276807840], [-651913, 4276677320], [-651943, 4276680904], [-651971, 4277725384], [-652079, 4277989064], [-652087, 4293590184], [-652217, 4294899944], [-532981, 2147516416], [-550369, 2155905024], [-642457, 3900743808], [-601231, 2863177728], [-641827, 3900735616]],
     [30887323656, 30069702976, 25811828928, 231962087968, 131097644288, 95834551808, 31139463360, 231936773184, 129260540544, 94627670032, 165356511376, 60431737104, 217231904768, 66035589216, 223338430992, 30737375760, 266556834688, 77645152544, 231970422928, 60151529520, 197770298370, 25791758353, 169114812488, 96746289664, 128871012384, 96169603088, 231950207552, 26446602760, 231970414752, 30069892104, 62411577360, 60171894856, 86979842304, 30887322640, 94824932368, 165658701828, 8589475839],
     16, 641827),
]


def combine(d1, d2, bound1, bound2):
    assert gcd(d1,d2) == 1
    k = d1*pow(d1,-1,d2)
    if 2*k > d1*d2:
        k -= d1*d2
    return abs(k)*bound1 + abs(1-k)*bound2


def main():
    routine = pinned(PATH)
    axes = list(map(abs,SIGNED_1+SIGNED_2))
    denominators, bounds = [], []
    fixed_rank = len(axes)
    for hidden,fibers,masks,scale,det in UPPER_DONORS:
        n = len(fibers)
        rows = [[i for i in range(n) if mask >> i & 1] for mask in masks]
        assert all(0 < mask < 1 << n for mask in masks)
        d = routine["donor"](hidden,fibers,rows,scale,det)
        denominators.append(d)
        bounds.append(Fraction(packet_bound(routine,fibers,rows,d),d))
        axes.extend(hidden)
        axes.extend(abs(p) for p,_ in fibers)
        fixed_rank += n+2
    assert len(set(axes)) == len(axes)
    du = denominators[0]
    upper = combine(*denominators,*bounds)
    upper, eu = blend(upper,bounds[0],du,0)
    lower_routine = pinned("scripts/ppi_lower_shield_pair.py")
    (dl,bl),(d2,b2) = [lower_routine["inspect"](c,m) for c,m in
                         [(SIGNED_1,MASKS_1),(SIGNED_2,MASKS_2)]]
    lower = combine(odd_part(dl),odd_part(d2),bl,b2)
    lower, el = blend(lower,bl,dl,max(dl&-dl,d2&-d2).bit_length()-1)
    E = odd_part(du*dl)
    dyadic = ((du&-du)*(dl&-dl)).bit_length()-1
    print("fixed axes",len(axes),"fixed rank bound",fixed_rank)
    select_sections(set(axes),E,dyadic,upper,lower,bounds[0],bl,eu,el,fixed_rank,2)


if __name__ == "__main__":
    main()
