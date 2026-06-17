from scipy.stats import kendalltau
import pandas as pd


def compute_kendall_tau(rankings):
    baseline = rankings["knowledge_baseline"]

    results = []

    for model, reps in rankings.items():
        if model == "knowledge_baseline":
            continue

        for rep_name, rep_ranks in reps.items():

            features = baseline.keys()

            baseline_vec = [baseline[f] for f in features]
            model_vec = [rep_ranks[f] for f in features]

            tau, pvalue = kendalltau(baseline_vec, model_vec)

            results.append(
                {
                    "model": model,
                    "rep": rep_name,
                    "kendall_tau": tau,
                    "pvalue": pvalue,
                }
            )

    return pd.DataFrame(results)




############## Robust learning from literature data ##################


"""
You are an expert in polymer physics, specifically conjugated polymer solution conformation and experimental factors governing Radius of Gyration (Rg).

Based on established literature, trained knowledge, and known physical principles, rank the following parameters by their influence on conjugated polymer solution conformation (Rg).

Task

Provide a strict ranking from 1 (most influential) to 12 (least influential).

Parameters to rank
Temperature SANS/SLS/DLS/SEC (K)
polymer dP
Mw (g/mol)
solvent dH
fp (polymer structural fingerprint / ECFP-based representation)
Concentration (mg/mL)
polymer dD
Xn (degree of polymerization)
solvent dP
PDI
polymer dH
solvent dD

Important definitions
fp = structural fingerprint representation of polymer (e.g., ECFP; encodes chemical structure)
Xn = degree of polymerization
dD, dP, dH = Hansen solubility parameters (dispersion, polarity, hydrogen bonding)
Temperature SANS/SLS/DLS/SEC (K): temperature of preparation
PDI: dispersity index
Mw: molecular weight
Concentration: Concentration of prepared sample in solution

Instructions
Use domain knowledge from polymer solution physics and scattering experiments (SANS/SLS/DLS/SEC).
Consider effects of:
chain expansion / contraction in solution
solvent quality
intermolecular interactions
concentration regimes (dilute/semi-dilute)
Do NOT leave any value as “None”.
Output format (STRICT)

Return only a valid JSON dictionary:

{
"Temperature SANS/SLS/DLS/SEC (K)": , 
"polymer dP": , 
"Mw (g/mol)": , 
"solvent dH": , 
"fp": , 
"Concentration (mg/ml)": ,
"polymer dD": ,
"Xn": ,
"solvent dP": ,
"PDI": ,
"polymer dH": ,
"solvent dD":
}

No explanation, no extra text.

"""
### this was done 3 time with opening new window


prompt1_robust_learning = {
  
  "knowledge_baseline": {
    "Xn": 1,
    "Mw (g/mol)": 2,
    "Concentration (mg/ml)": 3,
    "Temperature SANS/SLS/DLS/SEC (K)": 4,
    "polymer dD": 5,
    "solvent dD": 6,
    "polymer dP": 7,
    "solvent dP": 8,
    "polymer dH": 9,
    "solvent dH": 10,
    "fp": 11,
    "PDI": 12
  },

  "chatgpt55_extended_pro": {
    "rep_1": {
        "Temperature SANS/SLS/DLS/SEC (K)":11,
        "polymer dP":8,
        "Mw (g/mol)":1,
        "solvent dH":9,
        "fp":3,
        "Concentration (mg/ml)":4,
        "polymer dD":6,
        "Xn":2,
        "solvent dP":7,
        "PDI":12,
        "polymer dH":10,
        "solvent dD":5
    },
    "rep_2": {
        "Temperature SANS/SLS/DLS/SEC (K)": 12,
        "polymer dP": 9,
        "Mw (g/mol)": 2,
        "solvent dH": 7,
        "fp": 1,
        "Concentration (mg/ml)": 6,
        "polymer dD": 8,
        "Xn": 3,
        "solvent dP": 5,
        "PDI": 11,
        "polymer dH": 10,
        "solvent dD": 4
    },
    # "rep_3": {

    # }
  },

  "gemini31_pro_extended": {
    "rep_1": {
        "Temperature SANS/SLS/DLS/SEC (K)": 8,
        "polymer dP": 7,
        "Mw (g/mol)": 1,
        "solvent dH": 11,
        "fp": 2,
        "Concentration (mg/ml)": 9,
        "polymer dD": 5,
        "Xn": 3,
        "solvent dP": 6,
        "PDI": 10,
        "polymer dH": 12,
        "solvent dD": 4
    },
    "rep_2": {
        "Temperature SANS/SLS/DLS/SEC (K)": 7,
        "polymer dP": 10,
        "Mw (g/mol)": 1,
        "solvent dH": 11,
        "fp": 3,
        "Concentration (mg/ml)": 6,
        "polymer dD": 5,
        "Xn": 2,
        "solvent dP": 9,
        "PDI": 8,
        "polymer dH": 12,
        "solvent dD": 4
    },
    # "rep_3": {

    # }
  },

  "claude_opus48_max": {
    "rep_1": {
        "Temperature SANS/SLS/DLS/SEC (K)": 8,
        "polymer dP": 9,
        "Mw (g/mol)": 1,
        "solvent dH": 5,
        "fp": 3,
        "Concentration (mg/ml)": 6,
        "polymer dD": 10,
        "Xn": 2,
        "solvent dP": 4,
        "PDI": 12,
        "polymer dH": 11,
        "solvent dD": 7
    },
    "rep_2": {
        "Temperature SANS/SLS/DLS/SEC (K)": 4,
        "polymer dP": 10,
        "Mw (g/mol)": 1,
        "solvent dH": 7,
        "fp": 3,
        "Concentration (mg/ml)": 8,
        "polymer dD": 11,
        "Xn": 2,
        "solvent dP": 6,
        "PDI": 9,
        "polymer dH": 12,
        "solvent dD": 5
    },
    # "rep_3": {

    # }
  }
}


"""
According to your data that you're trained on and the papers about the Conjugated Polymer Solution Conformation and 
the experimental parameters that control Radius of Gyration (Rg), complete this ranking table:
PLS_Ranks = pd.DataFrame([{ 
"Temperature SANS/SLS/DLS/SEC (K)": None,
"polymer dP": None, 
"Mw (g/mol)": None,
"solvent dH": None,
"fp": None,
"Concentration (mg/ml)": None,
"polymer dD": None,
"Xn": None, 
"solvent dP": None, 
"PDI": None,
"polymer dH": None,
"solvent dD": None,
}], index=["expert_rank"])
from 1 (the most effective and predictive and controlling parameter) to least effective parameter.
FP can be consider as structure of polymer (fingerprint like ECFP encoded as single parameter,
Xn is degree of polymerization)

"""





prompt2_robust_learning = {
  "knowledge_baseline": {
    "Xn": 1,
    "Mw (g/mol)": 2,
    "Concentration (mg/ml)": 3,
    "Temperature SANS/SLS/DLS/SEC (K)": 4,
    "polymer dD": 5,
    "solvent dD": 6,
    "polymer dP": 7,
    "solvent dP": 8,
    "polymer dH": 9,
    "solvent dH": 10,
    "fp": 11,
    "PDI": 12
  },

  "chatgpt55_extended_pro": {
    "rep_1": {
      "Temperature SANS/SLS/DLS/SEC (K)": 4,
      "polymer dP": 10,
      "Mw (g/mol)": 2,
      "solvent dH": 11,
      "fp": 5,
      "Concentration (mg/ml)": 3,
      "polymer dD": 8,
      "Xn": 1,
      "solvent dP": 9,
      "PDI": 6,
      "polymer dH": 12,
      "solvent dD": 7
    },
    "rep_2": {
      "Temperature SANS/SLS/DLS/SEC (K)": 4,
      "polymer dP": 11,
      "Mw (g/mol)": 1,
      "solvent dH": 12,
      "fp": 6,
      "Concentration (mg/ml)": 2,
      "polymer dD": 9,
      "Xn": 3,
      "solvent dP": 7,
      "PDI": 5,
      "polymer dH": 10,
      "solvent dD": 8
    },
    "rep_3": {
      "Temperature SANS/SLS/DLS/SEC (K)": 4,
      "polymer dP": 9,
      "Mw (g/mol)": 2,
      "solvent dH": 11,
      "fp": 5,
      "Concentration (mg/ml)": 3,
      "polymer dD": 7,
      "Xn": 1,
      "solvent dP": 8,
      "PDI": 10,
      "polymer dH": 12,
      "solvent dD": 6
    }
  },

  "gemini31_pro_extended": {
    "rep_1": {
      "Temperature SANS/SLS/DLS/SEC (K)": 7,
      "polymer dP": 9,
      "Mw (g/mol)": 1,
      "solvent dH": 11,
      "fp": 3,
      "Concentration (mg/ml)": 6,
      "polymer dD": 5,
      "Xn": 2,
      "solvent dP": 8,
      "PDI": 10,
      "polymer dH": 12,
      "solvent dD": 4
    },
    "rep_2": {
      "Temperature SANS/SLS/DLS/SEC (K)": 6,
      "polymer dP": 8,
      "Mw (g/mol)": 1,
      "solvent dH": 9,
      "fp": 3,
      "Concentration (mg/ml)": 11,
      "polymer dD": 5,
      "Xn": 2,
      "solvent dP": 7,
      "PDI": 12,
      "polymer dH": 10,
      "solvent dD": 4
    },
    "rep_3": {
      "Mw (g/mol)": 1,
      "Xn": 2, 
      "fp": 3,
      "solvent dD": 4,
      "polymer dD": 5,
      "solvent dP": 6,
      "polymer dP": 7,
      "Temperature SANS/SLS/DLS/SEC (K)": 8,
      "Concentration (mg/ml)": 9,
      "PDI": 10,
      "solvent dH": 11,
      "polymer dH": 12
    }
  },

  "claude_opus48_max": {
    "rep_1": {
      "Temperature SANS/SLS/DLS/SEC (K)": 4,
      "polymer dP": 5,
      "Mw (g/mol)": 1,
      "solvent dH": 8,
      "fp": 9,
      "Concentration (mg/ml)": 3,
      "polymer dD": 11,
      "Xn": 2,
      "solvent dP": 6,
      "PDI": 12,
      "polymer dH": 7,
      "solvent dD": 10
    },
    "rep_2": {
      "Temperature SANS/SLS/DLS/SEC (K)": 10,
      "polymer dP": 8,
      "Mw (g/mol)": 1,
      "solvent dH": 6,
      "fp": 3,
      "Concentration (mg/ml)": 11,
      "polymer dD": 7,
      "Xn": 2,
      "solvent dP": 5,
      "PDI": 12,
      "polymer dH": 9,
      "solvent dD": 4
    },
    "rep_3": {
      "Temperature SANS/SLS/DLS/SEC (K)": 9,
      "polymer dP": 5,
      "Mw (g/mol)": 1,
      "solvent dH": 3,
      "fp": 7,
      "Concentration (mg/ml)": 8,
      "polymer dD": 11,
      "Xn": 2,
      "solvent dP": 4,
      "PDI": 12,
      "polymer dH": 6,
      "solvent dD": 10
    }
  }
}


### A prompet with paper included:
"""
According to your data that you're trained on and the papers about the Conjugated Polymer Solution Conformation and
the experimental parameters that control Radius of Gyration (Rg), complete this ranking table:
PLS_Ranks = pd.DataFrame([{ "Temperature SANS/SLS/DLS/SEC (K)": None,
"polymer dP": None, "Mw (g/mol)": None,
"solvent dH": None, "fp": None, "Concentration (mg/ml)": None,
"polymer dD": None, "Xn": None, "solvent dP": None, "PDI": None,
"polymer dH": None, "solvent dD": None, }], index=["expert_rank"])
from 1 (the most effective and predictive and controlling parameter)
to least effective parameter. FP can be consider as structure of polymer 
(fingerprint like ECFP encoded as single parameter) Xn is degree of polymerization.
PDI is polymer dispersity. Temperature SANS/SLS/DLS/SEC (K): preparation temperature. 
Mw: molecular weight dD, dH, dP: HSP parameters of solvent and polymer. 
You can look into papers like https://doi.org/10.1063/5.0303721 and similar!

"""

prompt3_robust_learning = {
  "knowledge_baseline": {
      "Xn": 1,
      "Mw (g/mol)": 2,
      "Concentration (mg/ml)": 3,
      "Temperature SANS/SLS/DLS/SEC (K)": 4,
      "polymer dD": 5,
      "solvent dD": 6,
      "polymer dP": 7,
      "solvent dP": 8,
      "polymer dH": 9,
      "solvent dH": 10,
      "fp": 11,
      "PDI": 12
  },

  "chatgpt55_extended_pro": {
    "rep_1": {
      "Temperature SANS/SLS/DLS/SEC (K)": 4,
      "polymer dP": 11,
      "Mw (g/mol)": 1,
      "solvent dH": 12,
      "fp": 8,
      "Concentration (mg/ml)": 2,
      "polymer dD": 9,
      "Xn": 3,
      "solvent dP": 6,
      "PDI": 5,
      "polymer dH": 10,
      "solvent dD": 7,
    },
    "rep_2": {
      "Temperature SANS/SLS/DLS/SEC (K)": 4,
      "polymer dP": 11,
      "Mw (g/mol)": 1,
      "solvent dH": 12,
      "fp": 6,
      "Concentration (mg/ml)": 2,
      "polymer dD": 9,
      "Xn": 3,
      "solvent dP": 7,
      "PDI": 5,
      "polymer dH": 10,
      "solvent dD": 8,
    },
    "rep_3": {
      "Temperature SANS/SLS/DLS/SEC (K)": 4,
      "polymer dP": 11,
      "Mw (g/mol)": 1,
      "solvent dH": 9,
      "fp": 6,
      "Concentration (mg/ml)": 2,
      "polymer dD": 10,
      "Xn": 3,
      "solvent dP": 8,
      "PDI": 5,
      "polymer dH": 12,
      "solvent dD": 7,
    }
  },

  "gemini31_pro_extended": {
    "rep_1": {
      "Mw (g/mol)": 1,
      "Xn": 2, 
      "fp": 3,
      "PDI": 4,
      "solvent dD": 5,
      "polymer dD": 6,
      "solvent dP": 7,
      "polymer dP": 8,
      "solvent dH": 9,
      "Concentration (mg/ml)": 10,
      "polymer dH": 11,
      "Temperature SANS/SLS/DLS/SEC (K)": 12
    },
    "rep_2": {
      "Mw (g/mol)": 1,
      "Xn": 2,
      "Concentration (mg/ml)": 3,
      "Temperature SANS/SLS/DLS/SEC (K)": 4,
      "solvent dD": 5,
      "solvent dP": 6,
      "polymer dD": 7,
      "polymer dP": 8,
      "solvent dH": 9,
      "polymer dH": 10,
      "PDI": 11,
      "fp": 12
    },
    "rep_3": {
      "Mw (g/mol)": 1,
      "Xn": 2, 
      "Concentration (mg/ml)": 3,
      "Temperature SANS/SLS/DLS/SEC (K)": 4,
      "solvent dH": 5, 
      "solvent dP": 6,
      "polymer dH": 7,
      "polymer dP": 8, 
      "solvent dD": 9, 
      "polymer dD": 10,
      "PDI": 11,
      "fp": 12
    }
  },

  "claude_opus48_max": {
    "rep_1": {
      "Temperature SANS/SLS/DLS/SEC (K)": 4,
      "polymer dP": 6,
      "Mw (g/mol)": 1,
      "solvent dH": 7,
      "fp": 11,
      "Concentration (mg/ml)": 3,
      "polymer dD": 10,
      "Xn": 2,
      "solvent dP": 5,
      "PDI": 12,
      "polymer dH": 8,
      "solvent dD": 9,
    },
    "rep_2": {
      "Temperature SANS/SLS/DLS/SEC (K)": 4,
      "polymer dP": 10,
      "Mw (g/mol)": 1,
      "solvent dH": 8,
      "fp": 5,
      "Concentration (mg/ml)": 3,
      "polymer dD": 9,
      "Xn": 2,
      "solvent dP": 7,
      "PDI": 12,
      "polymer dH": 11,
      "solvent dD": 6,
    },
    "rep_3": {
      "Temperature SANS/SLS/DLS/SEC (K)": 4,
      "polymer dP": 9,
      "Mw (g/mol)": 2,
      "solvent dH": 7,
      "fp": 5,
      "Concentration (mg/ml)": 3,
      "polymer dD": 12,
      "Xn": 1,
      "solvent dP": 8,
      "PDI": 11,
      "polymer dH": 10,
      "solvent dD": 6,
    }
  }
}

############## Beyond molecular structure_ critically assessing machine learning for designing organic photovoltaic ##################


"""
prompt 1

According to your data that you're trained on and the papers about the Conjugated Polymer Solution Conformation and
the experimental parameters that control Radius of Gyration (Rg), complete this ranking table:
PLS_Ranks = pd.DataFrame([{ "Temperature SANS/SLS/DLS/SEC (K)": None,
"polymer dP": None, "Mw (g/mol)": None,
"solvent dH": None, "fp": None, "Concentration (mg/ml)": None,
"polymer dD": None, "Xn": None, "solvent dP": None, "PDI": None,
"polymer dH": None, "solvent dD": None, }], index=["expert_rank"])
from 1 (the most effective and predictive and controlling parameter)
to least effective parameter. FP can be consider as structure of polymer 
(fingerprint like ECFP encoded as single parameter) Xn is degree of polymerization.
PDI is polymer dispersity. Temperature SANS/SLS/DLS/SEC (K): preparation temperature. 
Mw: molecular weight dD, dH, dP: HSP parameters of solvent and polymer. 
You can look into papers like https://doi.org/10.1063/5.0303721 and similar!


"""


prompt1_Beyond_molecular_structure = {
  "knowledge_baseline": {
      "HOMO_D (eV)": 1,
      "LUMO_A (eV)": 2,
      "Eg_D (eV)": 3,
      "Eg_A (eV)": 4,
      "Ehl_D (eV)": 5,
      "Ehl_A (eV)": 6,
      "HOMO_A (eV)": 7,
      "LUMO_D (eV)": 8,
      "fp_Donor": 9,
      "fp_Acceptor": 10,
      "D:A ratio (m/m)": 11,
      "temperature of thermal annealing": 12,
      "solvent additive conc. (% v/v)": 13,
      "HTL energy level (eV)": 14,
      "ETL energy level (eV)": 15
  },

  "chatgpt55_extended_pro": {
    "rep_1": {

    },
    "rep_2": {

    },
    "rep_3": {

    }
  },

  "gemini31_pro_extended": {
    "rep_1": {

    },
    "rep_2": {

    },
    "rep_3": {

    }
  },

  "claude_opus48_max": {
    "rep_1": {

    },
    "rep_2": {

    },
    "rep_3": {

    }
  }
}


prompt2_Beyond_molecular_structure = {
  "knowledge_baseline": {

  },

  "chatgpt55_extended_pro": {
    "rep_1": {

    },
    "rep_2": {

    },
    "rep_3": {

    }
  },

  "gemini31_pro_extended": {
    "rep_1": {

    },
    "rep_2": {

    },
    "rep_3": {

    }
  },

  "claude_opus48_max": {
    "rep_1": {

    },
    "rep_2": {

    },
    "rep_3": {

    }
  }
}

prompt3_Beyond_molecular_structure = {
  "knowledge_baseline": {

  },

  "chatgpt55_extended_pro": {
    "rep_1": {

    },
    "rep_2": {

    },
    "rep_3": {

    }
  },

  "gemini31_pro_extended": {
    "rep_1": {

    },
    "rep_2": {

    },
    "rep_3": {

    }
  },

  "claude_opus48_max": {
    "rep_1": {

    },
    "rep_2": {

    },
    "rep_3": {

    }
  }
}




if __name__ == "__main__":
    
    # prompt1_df = compute_kendall_tau(prompt1)
    # print(prompt1_df)

    # print("\n\n\n")
    # print("Kendall Tau Summary:")
    # print(prompt1_df.groupby("model")["kendall_tau"].agg(["mean", "std"]))

    # prompt2_df = compute_kendall_tau(prompt2)
    # print(prompt2_df)

    # print("\n\n\n")
    # print("Kendall Tau Summary:")
    # print(prompt2_df.groupby("model")["kendall_tau"].agg(["mean", "std"]))


    prompt3_df = compute_kendall_tau()
    print(prompt3_df)

    print("\n\n\n")
    print("Kendall Tau Summary:")
    print(prompt3_df.groupby("model")["kendall_tau"].agg(["mean", "std"]))


