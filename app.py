import streamlit as st
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw

# Load dataset
@st.cache_data
def load_data():
    df = pd.read_csv("PubChem_compound_text_antidepressant drugs.csv")
    df.dropna(subset=["Name", "SMILES"], inplace=True)
    return df

df = load_data()

# App title and sidebar
st.set_page_config(page_title="Antidepressant Drug Database", layout="wide")

# Tabs: Intro and Search
tab1, tab2 = st.tabs(["📘 Introduction", "🔍 Search Drugs"])

# Intro Screen
with tab1:
    st.title("🧠 Antidepressant Drug Database")
    st.markdown("""
    Welcome to the **Antidepressant Drug Database**, a specialized repository for chemical and pharmacological information on antidepressant compounds.  
    This app allows you to:
    - 🔍 Search for antidepressants by name  
    - 🧬 Visualize their 2D molecular structure using SMILES  
    - 📊 View detailed data for each compound  
    The data is derived from PubChem and curated into this compact format for academic and research use.
    
    ---
    """)

# Search & Display Tab
with tab2:
    st.header("🔍 Search Antidepressants")

    drug_names = sorted(df["Name"].unique())
    selected_name = st.selectbox("Select a drug by name", ["-- Select --"] + drug_names)

    if selected_name and selected_name != "-- Select --":
        selected_entry = df[df["Name"] == selected_name].iloc[0]
        st.subheader(f"📄 Details for: {selected_name}")

        # Show all details
        st.dataframe(selected_entry.to_frame().rename(columns={selected_entry.name: "Value"}))

        # Show molecule
        try:
            mol = Chem.MolFromSmiles(selected_entry["SMILES"])
            if mol:
                st.subheader("🧪 2D Molecular Structure")
                st.image(Draw.MolToImage(mol, size=(300, 300)))
            else:
                st.warning("Could not generate structure from SMILES.")
        except:
            st.error("Error rendering molecule. Invalid SMILES?")

    else:
        st.info("Please select a drug from the dropdown to view its details.")


