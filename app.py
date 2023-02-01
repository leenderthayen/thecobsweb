import streamlit as st
import csv
import numpy as np
import pandas as pd
import thecobs.spectral_functions as sf
import thecobs.constants as constants

from bokeh.models import DataRange1d
from bokeh.plotting import figure
from bokeh.palettes import Category10

apptitle = "The cob(s)web"

st.set_page_config(page_title=apptitle, page_icon=":eyeglasses:")

st.title("The Cob(s)web")

st.markdown("""
 Calculate beta spectrum shapes with state of the art methods
""")

st.sidebar.markdown("## Set the decay parameters")

select_event = st.sidebar.selectbox('How do you want to find data?',
                                    ['Manually', 'Query NNDC'])

@st.cache
def downloadAME():
        pass

if select_event == 'Manually':
    beta_type = st.sidebar.selectbox('Type of beta transition',
                                    ['Beta-', 'Beta+'])

    sl_z = st.sidebar.number_input('Proton number', min_value=1, max_value=120, step=1, help="Set the proton number of the initial state")
    sl_a = st.sidebar.number_input('Mass number', min_value=sl_z, max_value=120, step=1, help="Set the mass number of the initial state")
    sl_r = st.sidebar.number_input('Radius', min_value=0.01, max_value=100., step=0.1, value=1.2, help="Set the rms nuclear radius in fm")
    sl_e0 = st.sidebar.number_input('Endpoint energy', min_value=0., max_value=20., value=1., step=0.001, help="Set the endpoint energy in MeV")

    z = int(sl_z)
    a = int(sl_a)
    e0 = float(sl_e0)
    r = float(sl_r)*(5./3.)**0.5*1e-15/constants.NATURALLENGTH

else:
    str_iso = st.sidebar.text_input('Isotope', placeholder='Ex. 6He')

    z = 2
    a = 6
    e0 = 3.5
    e_step = 0.001
    r = 0.01

sl_e_step = st.sidebar.number_input('Energy step', min_value=0.001, max_value=1., value=0.001, step=0.001)
e_step = float(sl_e_step)

@st.cache
def calculateSpectrum(Z, A, R, E0, E_step):
    E = np.arange(0, E0, E_step)

    W0 = 1 + E0/constants.ELECTRON_MASS_C2

    W = 1 + E/constants.ELECTRON_MASS_C2

    ph = sf.phase_space(W, W0)
    f = sf.fermi_function(Z, W, R)
    l0 = sf.finite_size_L0(Z, W, R)
    u = sf.finite_size_U_fermi(Z, W)
    rc = sf.radiative_correction(Z, W, W0, R)

    sp = ph*f*l0*u*rc

    comb = np.stack((E, W, sp, ph, f, l0, u, rc), axis=1)

    df = pd.DataFrame(comb, columns = ['Energy', 'W', 'Spectrum', 'PhaseSpace', 'FermiFunction', 'L0', 'U', 'RadiativeCorrections'])

    return df

@st.cache
def convert_df(df):
    # IMPORTANT: Cache the conversion to prevent computation on every rerun
    return df.to_csv(index=False).encode('utf-8')


df = calculateSpectrum(z, a, r, e0, e_step)

st.subheader('Electron spectrum')

p = figure(title="Electron spectrum", x_axis_label='Kinetic energy [MeV]', y_axis_label='dN/dE')

p.line(df['Energy'], df['Spectrum'], legend_label='All corrections', line_width=2, muted_alpha=0.2, color=Category10[3][0])
p.line(df['Energy'], df['PhaseSpace'], legend_label='Phase Space', line_width=2, muted_alpha=0.2, color=Category10[3][1])

p.legend.click_policy = "mute"

st.bokeh_chart(p, use_container_width=True)

csv = convert_df(df)

st.download_button(
    label="Download data as CSV",
    data=csv,
    file_name='spectrum.csv',
    mime='text/csv',
)


# -- Notes on whitening
with st.expander("More details"):
    st.markdown("""
 You are looking at the probability distribution for the kinetic energy of the beta particle after (nuclear) beta decay.
 The different corrections can be plotted individually by clicking the elements in the legend
""")
    if st.checkbox('Show individual corrections'):
        pcorr = figure(title='Corrections', x_axis_label='Kinetic energy [MeV]', y_axis_label='Correction', y_range=DataRange1d(only_visible=True))
        i = 0
        dfCorr = df[['FermiFunction', 'L0', 'U', 'RadiativeCorrections']]
        for column in dfCorr:
            pcorr.line(df['Energy'], df[column], legend_label=column, color=Category10[len(dfCorr.columns)][i], line_width=2)
            i+=1
        pcorr.legend.click_policy = "hide"

        st.bokeh_chart(pcorr, use_container_width=True)

    if st.checkbox('Show raw data'):
        df

