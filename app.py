import streamlit as st
import pandas as pd
import requests
from st_aggrid import AgGrid, GridOptionsBuilder

@st.cache_data
def load_data():
    return pd.read_csv("PubChem_compound_text_antidepressant drugs.csv").fillna("N/A")

df = load_data()

# Set page config
st.set_page_config(page_title="Antidepressant Drug Info", layout="wide")
st.title("Antidepressant Drug Information Portal")

# Sidebar for search and filters
with st.sidebar:
    st.header("Search Filters")
    query = st.text_input("Enter Drug Name or Synonym")
    selected_class = st.selectbox("Filter by Drug Class", 
                                options=["All"] + sorted(df["Drug_Class"].unique().tolist()))
    
    search_button = st.button("Search Drugs")

# Main content area
if search_button:
    # Apply filters
    filtered = df[df["Name"].str.contains(query, case=False, na=False) |
                 df["Synonyms"].str.contains(query, case=False, na=False)]
    
    if selected_class != "All":
        filtered = filtered[filtered["Drug_Class"] == selected_class]

    if filtered.empty:
        st.warning("No results found")
    else:
        # Create PubChem links
        filtered['PubChem_Link'] = filtered['PubChem_CID'].apply(
            lambda x: f'https://pubchem.ncbi.nlm.nih.gov/compound/{x}'
        )
        
        # Configure AgGrid
        gb = GridOptionsBuilder.from_dataframe(
            filtered[["Name", "Molecular_Formula", "Molecular_Weight", "Drug_Class", "PubChem_Link"]]
        )
        gb.configure_column("PubChem_Link", 
                           headerName="PubChem Link",
                           cellRenderer=lambda x: f'<a href="{x}" target="_blank">🔗 View</a>',
                           width=120)
        gb.configure_selection('single', use_checkbox=False)
        grid_options = gb.build()
        
        # Display results grid
        st.header("Search Results")
        grid_response = AgGrid(
            filtered,
            gridOptions=grid_options,
            height=300,
            width='100%',
            theme='streamlit',
            update_mode='SELECTION_CHANGED',
            fit_columns_on_grid_load=True
        )
        
        # Show detailed view when a row is selected
        if grid_response['selected_rows']:
            selected_drug = grid_response['selected_rows'][0]
            st.divider()
            
            # Create tabs for different information categories
            tab1, tab2, tab3, tab4 = st.tabs(["Basic Info", "Chemical Properties", 
                                              "Safety & Effects", "Vendor & References"])
            
            with tab1:
                col1, col2 = st.columns([1, 1])
                with col1:
                    st.subheader("Identification")
                    st.markdown(f"**Name:** {selected_drug.get('Name', 'N/A')}")
                    st.markdown(f"**Synonyms:** {selected_drug.get('Synonyms', 'N/A')}")
                    st.markdown(f"**Drug Class:** {selected_drug.get('Drug_Class', 'N/A')}")
                    
                with col2:
                    st.subheader("Molecular Properties")
                    st.markdown(f"**Formula:** {selected_drug.get('Molecular_Formula', 'N/A')}")
                    st.markdown(f"**Weight:** {selected_drug.get('Molecular_Weight', 'N/A')}")
                    st.markdown(f"**XLogP:** {selected_drug.get('XLogP', 'N/A')}")
            
            with tab2:
                col1, col2 = st.columns([1, 1])
                with col1:
                    st.subheader("Structural Information")
                    st.markdown(f"**IUPAC Name:** {selected_drug.get('IUPAC_Name', 'N/A')}")
                    st.markdown(f"**Canonical SMILES:** `{selected_drug.get('Canonical_SMILES', 'N/A')}`")
                    st.markdown(f"**InChI Key:** `{selected_drug.get('InChIKey', 'N/A')}`")
                    
                with col2:
                    st.subheader("Physicochemical Properties")
                    st.markdown(f"**Hydrogen Bond Donor:** {selected_drug.get('Hydrogen_Bond_Donor_Count', 'N/A')}")
                    st.markdown(f"**Hydrogen Bond Acceptor:** {selected_drug.get('Hydrogen_Bond_Acceptor_Count', 'N/A')}")
                    st.markdown(f"**Rotatable Bonds:** {selected_drug.get('Rotatable_Bond_Count', 'N/A')}")
            
            with tab3:
                col1, col2 = st.columns([1, 1])
                with col1:
                    st.subheader("Clinical Information")
                    st.markdown(f"**Indication:** {selected_drug.get('Indication', 'N/A')}")
                    st.markdown(f"**Mechanism of Action:** {selected_drug.get('Mechanism_of_Action', 'N/A')}")
                    
                with col2:
                    st.subheader("Safety Profile")
                    st.markdown(f"**Side Effects:** {selected_drug.get('Side_Effects', 'N/A')}")
                    st.markdown(f"**Toxicity:** {selected_drug.get('Toxicity', 'N/A')}")
                    st.markdown(f"**Precautions:** {selected_drug.get('Precautions', 'N/A')}")
            
            with tab4:
                col1, col2 = st.columns([1, 1])
                with col1:
                    st.subheader("Vendor Information")
                    st.markdown(f"**Supplier:** {selected_drug.get('Vendor', 'N/A')}")
                    st.markdown(f"**Catalog Number:** {selected_drug.get('Catalog_Number', 'N/A')}")
                    
                with col2:
                    st.subheader("References")
                    st.markdown(f"**PubChem:** [View Details]({selected_drug.get('PubChem_Link', '')})")
                    st.markdown(f"**ChEMBL ID:** {selected_drug.get('ChEMBL_ID', 'N/A')}")
                    st.markdown(f"**Therapeutic Targets DB:** {selected_drug.get('TTD_ID', 'N/A')}")

# Initial state before search
else:
    st.info("Enter search criteria and click 'Search Drugs' to begin")
    with st.expander("Dataset Preview"):
        st.dataframe(df.head(10))
