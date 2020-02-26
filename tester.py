import pandas as pd

df = pd.read_csv("Svarog data")

def simple_cal(track, formula_data):
    print(track)
    print(formula_data)
    return len(track), len(formula_data)