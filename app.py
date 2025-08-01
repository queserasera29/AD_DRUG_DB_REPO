import streamlit as st
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
from io import BytesIO

# Load data
@st.cache_data
def load_data():
    df = pd.read_csv("PubChem_compound_text_antidepressant drugs.csv")
    df = df.dropna(subset=["Name", "SMILES"])
    return df

df = load_data()

# Function to render molecule from SMILES
def render_mol(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        img = Draw.MolToImage(mol, size=(300, 300))
        buf = BytesIO()
        img.save(buf, format="PNG")
        return buf.getvalue()
    return None

# Streamlit app layout
st.set_page_config(page_title="Antidepressant Drug Database", layout="wide")

# Sidebar navigation
page = st.sidebar.radio("Navigation", ["🏠 Home", "🔍 Search Drug"])

# Page: Intro
if page == "🏠 Home":
    st.title("💊 Antidepressant Drug Database")
    st.markdown(
        """
        Welcome to the **Antidepressant Drug Database** — a curated collection of compounds commonly used 
        in the treatment of depression and related mental health disorders. This tool allows you to:
        - Search drugs by name  
        - View molecular structures  
        - Explore pharmacological and chemical properties  
        
        This database is built using data from **PubChem**, tailored specifically for antidepressants.
        """
    )
    st.image("https://images.unsplash.com/photo-1580281657527-47d4f7e2b3d6", caption="For academic & research use only", use_column_width=True)

# Page: Drug Search
elif page == "🔍 Search Drug":
    st.title("🔍 Search Antidepressant Drug")

    query = st.text_input("Enter Drug Name").strip().lower()

    if query:
        results = df[df["Name"].str.lower().str.contains(query)]
        if not results.empty:
            for _, row in results.iterrows():
                st.subheader(row["Name"])
                
                col1, col2 = st.columns([1, 2])
                with col1:
                    mol_img = render_mol(row["SMILES"])
                    if mol_img:
                        st.image(mol_img, caption="2D Structure", use_column_width=False)
                    else:
                        st.error("Could not render structure.")
                with col2:
                    st.markdown("**Details:**")
                    for col in df.columns:
                        st.markdown(f"**{col}**: {row[col]}")
                st.markdown("---")
        else:
            st.warning("No results found. Please check your spelling.")
    else:
        st.info("Enter a drug name in the search box above.")

