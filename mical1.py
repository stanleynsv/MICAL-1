
import mols2grid
import pandas as pd
import streamlit as st
import streamlit.components.v1 as components
from rdkit import Chem
from rdkit.Chem.Descriptors import ExactMolWt, MolLogP

st.title("MICAL-1 Inhibitors")


@st.cache_data()
def download_dataset():
    """Loads once then cached for subsequent runs"""
    df = pd.read_excel(
        "MICAL1_inhibitors.xlsx"
    ).dropna()
    return df

# Calculate descriptors
def calc_mw(smiles_string):
    """Given a smiles string (ex. C1CCCCC1), calculate and return the molecular weight"""
    mol = Chem.MolFromSmiles(smiles_string)
    return ExactMolWt(mol)

def calc_logp(smiles_string):
    """Given a smiles string (ex. C1CCCCC1), calculate and return the LogP"""
    mol = Chem.MolFromSmiles(smiles_string)
    return MolLogP(mol)

# Copy the dataset so any changes are not applied to the original cached version
df = download_dataset().copy()
df["MW"] = df.apply(lambda x: calc_mw(x["SMILES"]), axis=1)
df["LogP"] = df.apply(lambda x: calc_logp(x["SMILES"]), axis=1)

# Sidebar panel
st.sidebar.header('Set parameters')
st.sidebar.write('*Note: Display compounds having values less than the following thresholds*')

# weight_cutoff = st.sidebar.slider("Molecular weight", 0, 1000, (300, 600), step=10)
inhibitor_efficiency_half_cutoff = st.sidebar.slider(
    label="0.5 µM inhibitor efficiency (%)",
    min_value=0.0,
    max_value=200.0,
    value=100.0,
    step=0.1,
)
inhibitor_efficiency_5_cutoff = st.sidebar.slider(
    label="5 µM inhibitor efficiency (%)",
    min_value=0.0,
    max_value=200.0,
    value=100.0,
    step=0.1,
)
inhibitor_efficiency_50_cutoff = st.sidebar.slider(
    label="50 µM inhibitor efficiency (%)",
    min_value=0.0,
    max_value=200.0,
    value=100.0,
    step=0.1,
)
# weight_cutoff = st.sidebar.slider(
#     label="Molecular weight",
#     min_value=0,
#     max_value=1000,
#     value=500,
#     step=10,
# )
# logp_cutoff = st.sidebar.slider(
#     label="LogP",
#     min_value=-10,
#     max_value=10,
#     value=5,
#     step=1,
# )

df_result = df[df["0.5 µM inhibitor efficiency (%)"] < inhibitor_efficiency_half_cutoff]
df_result2 = df_result[df_result["5 µM inhibitor efficiency (%)"] < inhibitor_efficiency_5_cutoff]
df_result3 = df_result2[df_result2["50 µM inhibitor efficiency (%)"] < inhibitor_efficiency_50_cutoff]

# df_result2 = df_result[df_result["MW"] < weight_cutoff]
# df_result3 = df_result2[df_result2["LogP"] < logp_cutoff]

st.write(df_result3.shape)
st.write(df_result3)

raw_html = mols2grid.display(df_result3,
                             subset=["Number", "img", "0.5 µM inhibitor efficiency (%)", "5 µM inhibitor efficiency (%)", "50 µM inhibitor efficiency (%)"],
                             mapping={"SMILES": "Smiles", "Number": "Number", "0.5 µM inhibitor efficiency (%)": "0.5 µM inhibitor efficiency (%)"})._repr_html_()
components.html(raw_html, width=900, height=1200, scrolling=False)


