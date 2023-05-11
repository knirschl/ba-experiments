import os
import sys

def get_model(subst_model):
  return subst_model.split("+")[0]

def get_gamma_rates(subst_model):
  if ("G" in subst_model.split("+")):
    return 4
  else:
    return 1

def is_dna(subst_model):
  dna_models = ["JC", "GTR"]
  return get_model(subst_model) in dna_models

def get_raxml_model(subst_model):
  subst_model = subst_model.replace("+I", "+IC") # empirical invariant sites
  res = get_model(subst_model)
  if (res == "POISSON"):
    res = "PROTGTR{"
    res += "/".join(["1"] * 190)
    res += "}"
    res += "+FE"
  rest = subst_model.split("+")[1:]
  for r in rest:
    res += "+" + r
  return res