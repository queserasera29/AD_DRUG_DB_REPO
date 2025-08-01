import streamlit as st
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw

# Load the data
df = pd.read_csv("PubChem_compound_text_antidepressant drugs.csv")

# Sidebar navigation
page = st.sidebar.selectbox("Go to", ["Home", "Drug Search"])

# Home Page
if page == "Home":
    st.title("🧠 Antidepressant Drug Database")
    st.markdown("""
    Welcome to the **Antidepressant Drug Knowledgebase**!

    This is a mini PubChem-style resource focused exclusively on **antidepressants**.

    ### 🔍 Features:
    - Search drugs by name
    - View detailed compound information
    - See 2D molecular structures rendered from SMILES

    **Use the sidebar** to begin searching.
    """)

# Drug Search Page
elif page == "Drug Search":
    st.title("🔬 Explore Antidepressant Compounds")

    search_query = st.text_input("Search by Drug Name")

    # Filter data
    if search_query:
        filtered_df = df[df["Name"].str.contains(search_query, case=False, na=False)]
    else:
        filtered_df = df.copy()

    # Show results
    if not filtered_df.empty:
        for idx, row in filtered_df.iterrows():
            st.subheader(row.get("Name", "Unknown Name"))

            # Drug details
            st.markdown(f"""
            **Formula:** {row.get("Formula", "N/A")}  
            **Molecular Weight:** {row.get("Molecular_Weight", "N/A")}  
            **XLogP:** {row.get("XLogP", "N/A")}  
            **Synonyms:** {row.get("Synonyms", "N/A")}
            """)

            # 2D structure rendering
            smiles = row.get("SMILES", "")
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                img = Draw.MolToImage(mol, size=(300, 300))
                st.image(img, caption="2D Structure", use_column_width=False)
            else:
                st.warning("No valid SMILES available for 2D structure.")
            
            st.markdown("---")
    else:
        st.info("No matching drugs found. Try a different search term.")
