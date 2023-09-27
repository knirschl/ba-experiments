setup_10 = {
    'SPECIES': {
        '15': {'BA+FM': 60.035879153378154, 'FastME': 0.7870688779013497, 'RAxML': 27.406733741565628, 'GeneRax': 34.18040821017051}, '100': {'BA+FM': 133.7752560948839, 'FastME': 42.739728236684996, 'RAxML': 469.8387166213016, 'GeneRax': 616.557478272185}, '50': {'BA+FM': 75.47978641427294, 'FastME': 4.952131560870579, 'RAxML': 155.43138814459041, 'GeneRax': 190.41438170355192}, '75': {'BA+FM': 108.21290671886229, 'FastME': 31.06383581794038, 'RAxML': 289.4818441892157, 'GeneRax': 367.5987963652124}, '25': {'BA+FM': 65.22204767867011, 'FastME': 1.7916462372760384, 'RAxML': 53.94659502895511, 'GeneRax': 68.17047372886113}},
    'SITES': {
        '50': {'BA+FM': 66.13684868082709, 'FastME': 4.618998634571931, 'RAxML': 41.65644736435949, 'GeneRax': 58.45018747631384}, '500': {'BA+FM': 62.13679006634926, 'FastME': 1.9516922308474172, 'GeneRax': 119.94970300489543, 'RAxML': 142.94258653631016}, '250': {'BA+FM': 67.74444105734631, 'FastME': 1.6398659214681508, 'GeneRax': 87.13769552415731, 'RAxML': 87.35209701012592}, '100': {'BA+FM': 65.22204767867011, 'FastME': 1.7916462372760384, 'RAxML': 53.94659502895511, 'GeneRax': 68.17047372886113}},
    'BRALEN': {
        '0.01': {'BA+FM': 52.72354474177166, 'RAxML': 22.826985955238342, 'GeneRax': 67.80520010724359, 'FastME': 92.56689747255676}, '10.0': {'BA+FM': 62.93338037875233, 'FastME': 2.08930043298371, 'RAxML': 69.07263780126766, 'GeneRax': 120.68536375979988}, '100.0': {'BA+FM': 61.6522389601688, 'FastME': 2.0602877359001006, 'RAxML': 93.27139320178908, 'GeneRax': 194.07654236774056}, '0.1': {'BA+FM': 57.94453232385675, 'FastME': 3.363065534708451, 'RAxML': 35.249844808967744, 'GeneRax': 60.128633083129415}, '1.0': {'BA+FM': 65.22204767867011, 'FastME': 1.7916462372760384, 'RAxML': 53.94659502895511, 'GeneRax': 68.17047372886113}},
    'DUPLOS': {
        '0.5': {'BA+FM': 49.15823682960199, 'FastME': 1.6422387994065577, 'RAxML': 38.59385104325353, 'GeneRax': 72.82850963972052}, '0.0': {'BA+FM': 53.55077088244107, 'FastME': 1.4719196095758555, 'RAxML': 25.725944526341497, 'GeneRax': 59.30984746193399}, '2.0': {'BA+FM': 66.4727047535838, 'FastME': 14.987837861995308, 'RAxML': 109.02087859961452, 'GeneRax': 151.84681175436293}, '1.0': {'BA+FM': 65.22204767867011, 'FastME': 1.7916462372760384, 'RAxML': 53.94659502895511, 'GeneRax': 68.17047372886113}, '3.0': {'BA+FM': 116.35848952373678, 'FastME': 143.88378606629126, 'RAxML': 179.27731640728152, 'GeneRax': 217.06478972287522}}}   
setup_80 = {
    'SPECIES': {
        '15': {'BA+FM': 574.4330004554004}, '100': {'BA+FM': 626.5216008284277}, '50': {'BA+FM': 670.9986095043145}, '75': {'BA+FM': 1114.3331069737571}, '25': {'BA+FM': 567.5231903137181}},
    'SITES': {
        '50': {'BA+FM': 689.9190535593515}, '500': {'BA+FM': 588.0722097636071}, '250': {'BA+FM': 603.4600196775765}, '100': {'BA+FM': 567.5231903137181}},
    'BRALEN': {
        '0.01': {'BA+FM': 749.2525479183454}, '10.0': {'BA+FM': 592.5197143912316}, '100.0': {'BA+FM': 591.9491177304586}, '0.1': {'BA+FM': 530.0506950901813}},
    'DUPLOS': {
        '0.5': {'BA+FM': 528.2639875460151}, '0.0': {'BA+FM': 370.15611397216617}, '2.0': {'BA+FM': 762.1761020249392}, '1.0': {'BA+FM': 567.5231903137181}, '3.0': {'BA+FM': 1098.4811325394344}}}
# setup_ = {dataset : {variant : {method : runtime}}}


BASE = setup_10["DUPLOS"]["1.0"]
bafm80 = [setup_80[d][v]["BA+FM"] for d in setup_80 for v in setup_80[d]]
print("AVG-RT-BAFM80/3", sum(bafm80) / len(bafm80))
all10 = [setup_10[d][v]["BA+FM"] for d in setup_10 for v in setup_10[d]]
print("AVG-RT-ALL10", sum(all10) / len(all10))
fm = [setup_10[d][v]["FastME"] for d in setup_10 for v in setup_10[d]]
print("AVG-RT-FM", sum(fm) / len(fm))
rm = [setup_10[d][v]["RAxML"] for d in setup_10 for v in setup_10[d]]
print("AVG-RT-RM", sum(rm) / len(rm))
gr = [setup_10[d][v]["GeneRax"] for d in setup_10 for v in setup_10[d]]
print("AVG-RT-GR", sum(gr) / len(gr))
print("REL-DIFF-FACTOR-AVG-RT-BAFM80/3-ALL10", (sum(bafm80) / len(bafm80)) / (sum(all10) / len(all10)), ";inv:", (sum(all10) / len(all10)) / (sum(bafm80) / len(bafm80)))
print("REL-TIME-DIFF-ALL-FM", (sum(all10) / len(all10)) / (sum(fm) / len(fm)), ";inv:", (sum(fm) / len(fm)) / (sum(all10) / len(all10)))
all10 = [setup_10[d][v]["BA+FM"] for d in ["BRALEN", "SITES"] for v in setup_10[d]]
print("AVG-RT-ALL10-BL/SITES", sum(all10) / len(all10))
fm = [setup_10[d][v]["FastME"] for d in ["BRALEN", "SITES"] for v in setup_10[d]]
print("AVG-RT-FM-BL/SITES", sum(fm) / len(fm))
rm = [setup_10[d][v]["RAxML"] for d in ["BRALEN", "SITES"] for v in setup_10[d]]
print("AVG-RT-RM-BL/SITES", sum(rm) / len(rm))
gr = [setup_10[d][v]["GeneRax"] for d in ["BRALEN", "SITES"] for v in setup_10[d]]
print("AVG-RT-GR-BL/SITES", sum(gr) / len(gr))
print("RT-FM-BL0.01", setup_10["BRALEN"]["0.01"]["FastME"])
print("RT-FM-BASE", BASE["FastME"])
print("RT-ALL-BASE", BASE["BA+FM"])
print("REL-RT-DIFF-FM-BASE/SPECIES100", BASE["FastME"] / setup_10["SPECIES"]["100"]["FastME"], ";inv:", setup_10["SPECIES"]["100"]["FastME"] / BASE["FastME"])
print("RT-FM-SPECIES100", setup_10["SPECIES"]["100"]["FastME"])
print("REL-RT-DIFF-ALL-BASE/SPECIES100", BASE["BA+FM"] / setup_10["SPECIES"]["100"]["BA+FM"], ";inv:", setup_10["SPECIES"]["100"]["BA+FM"] / BASE["BA+FM"])
print("REL-RT-DIFF-GR-BASE/SPECIES100", BASE["GeneRax"] / setup_10["SPECIES"]["100"]["GeneRax"], ";inv:", setup_10["SPECIES"]["100"]["GeneRax"] / BASE["GeneRax"])
print("REL-RT-DIFF-RM-BASE/SPECIES100", BASE["RAxML"] / setup_10["SPECIES"]["100"]["RAxML"], ";inv:", setup_10["SPECIES"]["100"]["RAxML"] / BASE["RAxML"])
print("RT-ALL-SPECIES100", setup_10["SPECIES"]["100"]["BA+FM"])
print("RT-GR-DUPLOS3", setup_10["DUPLOS"]["3.0"]["GeneRax"])
rm = [setup_10[d][v]["RAxML"] for d in ["DUPLOS"] for v in setup_10[d]]
gr = [setup_10[d][v]["GeneRax"] for d in ["DUPLOS"] for v in setup_10[d]]
all10 = [setup_10[d][v]["BA+FM"] for d in ["DUPLOS"] for v in setup_10[d]]
print("AVG-RT-DIFF-AVGGR/RM-ALL-DUPLOS", (sum(rm+gr) / len(rm+gr)) / (sum(all10) / len(all10)))
print("ABS-TIME-DIFF-ALL-GR-SITES50", setup_10["SITES"]["50"]["BA+FM"] - setup_10["SITES"]["50"]["GeneRax"])
print("ABS-TIME-DIFF-ALL-GR-BL100", setup_10["BRALEN"]["100.0"]["BA+FM"] - setup_10["BRALEN"]["100.0"]["GeneRax"])
print("ABS-TIME-DIFF-ALL-RM-BL100", setup_10["BRALEN"]["100.0"]["BA+FM"] - setup_10["BRALEN"]["100.0"]["RAxML"])
print("ABS-TIME-DIFF-ALL-GR-SPECIES100", setup_10["SPECIES"]["100"]["BA+FM"] - setup_10["SPECIES"]["100"]["GeneRax"])
print("ABS-TIME-DIFF-ALL-RM-SPECIES100", setup_10["SPECIES"]["100"]["BA+FM"] - setup_10["SPECIES"]["100"]["RAxML"])
print(setup_10["SPECIES"]["100"]["BA+FM"])
print(setup_10["SPECIES"]["100"]["GeneRax"])