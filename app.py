import streamlit as st
import csv
import numpy as np
import pandas as pd
import re
import thecobs.SpectralFunctions as sf
from thecobs.Constants import *
import thecobs.Functions as fu
import thecobs.Screening

from bokeh.models import DataRange1d
from bokeh.plotting import figure
from bokeh.palettes import Category10

from io import StringIO
import requests

apptitle = "The cob(s)web"

st.set_page_config(page_title=apptitle, page_icon=":eyeglasses:")

st.title("The Cob(s)web")

st.markdown("""
 Calculate beta spectrum shapes with state of the art methods
""")

st.sidebar.markdown("## Set the decay parameters")

select_event = st.sidebar.selectbox('How do you want to find data?',
                                    ['Database', 'Manually'])

def approxNuclRadius(Z, A, Elton=False):
    N = A-Z
    R = 0
    if Elton:
        R = max(1, fu.getEltonNuclearRadius(A))
    else:
        rA = 1.282
        b = 0.342
        R = rA*(1-b*(N-Z)/A)*A**(1/3)
    return R*(3/5)**0.5

@st.cache
def downloadAME():
    names = ['N', 'Z', 'A', 'name', 'mass']
    data = []
    url = 'https://www-nds.iaea.org/amdc/ame2020/mass_1.mas20.txt'
    r = requests.get(url).content
    f = StringIO(r.decode('utf-8'))
    skipLines = 36
    currentLine = 0
    for line in f:
        currentLine += 1
        if currentLine <= skipLines:
            continue
        else:
            nMinusz = int(line[1:4])
            n = int(line[4:9])
            z = int(line[9:14])
            a = int(line[14:19])
            name = line[19:23]
            atomicMassString = line[106:121]
            atomicMassString = atomicMassString[:3] + atomicMassString[4:]
            atomicMassString = atomicMassString.replace('#', '.')
            atomicMass = float(atomicMassString)*1e-6*AMU_MASS_KEV
            data.append([n, z, a, name, atomicMass])
    return pd.DataFrame(data, columns=names)

decay_type = st.sidebar.selectbox('Type of beta transition',
                                    ['Beta-', 'Beta+'])

beta_type = st.sidebar.selectbox('Type of allowed transition',
                                    ['Gamow-Teller', 'Fermi'])

Qvalue = 1000.

if select_event == 'Manually':
    sl_z = st.sidebar.number_input('Proton number', min_value=1, max_value=120, step=1, help="Set the proton number of the initial state")
    sl_a = st.sidebar.number_input('Mass number', min_value=sl_z, max_value=120, step=1, value=int(fu.stableA(sl_z)), help="Set the mass number of the initial state")
    sl_r = st.sidebar.number_input('Radius', min_value=0.01, max_value=100., step=0.1, value=approxNuclRadius(sl_z, sl_a), help="Set the rms nuclear radius in fm")

    z = int(sl_z)
    a = int(sl_a)
    r = float(sl_r)*(5/3)**0.5*1e-15/NATURAL_LENGTH

else:
    dfAME = downloadAME()

    str_iso = st.sidebar.text_input('Isotope', value='1N')

    m = re.search(r'\d+', str_iso)
    if m:
        a = int(m.group(0))
        z = atoms.index(str(str_iso).replace(m.group(0), '').strip())
    else:
        st.error('Not a valid isotope name. Expecting something like 6He or 45Ca.')

    sl_r = st.sidebar.number_input('Radius', min_value=0.01, max_value=100., step=0.1, value=approxNuclRadius(z, a), help="Set the rms nuclear radius in fm")

    r = sl_r*(5./3.)**0.5*1e-15/NATURAL_LENGTH

    if decay_type == 'Beta-':
        Qvalue = dfAME.loc[(dfAME['Z'] == z) & (dfAME['A'] == a), 'mass'].values[0]-dfAME.loc[(dfAME['Z'] == (z+1)) & (dfAME['A'] == a), 'mass'].values[0]
    else:
        Qvalue = dfAME.loc[(dfAME['Z'] == z) & (dfAME['A'] == a), 'mass'].values[0]-dfAME.loc[(dfAME['Z'] == (z-1)) & (dfAME['A'] == a), 'mass'].values[0]-2*ELECTRON_MASS_KEV

sl_e0 = st.sidebar.number_input('Endpoint energy', min_value=0., max_value=20e3, value=Qvalue, step=1., help="Set the endpoint energy in keV")
sl_e_step = st.sidebar.number_input('Energy step', min_value=0.1, max_value=1000., value=1., step=1., help="Set the step energy in keV")

@st.cache
def calculateSpectrum(Z, A, R, E0, E_step, beta_type):
    E = np.arange(0.1, E0, E_step)

    W0 = 1 + E0/ELECTRON_MASS_KEV

    W = 1 + E/ELECTRON_MASS_KEV

    ph = sf.phase_space(W, W0)
    f = sf.fermi_function(W, Z, R)
    l0 = sf.finite_size_L0(W, Z, R)
    u = sf.finite_size_U_fermi(W, Z)
    rc = sf.radiative_correction(W, Z, W0, R)

    if beta_type == 'Fermi':
        rec = sf.recoil_fermi(W, W0, A)
        recCoul = sf.recoil_Coulomb_fermi(W, Z, W0, A)
        C = sf.shape_factor_fermi(W, Z, W0, R)
    elif beta_type == 'Gamow-Teller':
        rec = sf.recoil_gamow_teller(W, W0, A)
        recCoul = sf.recoil_Coulomb_gamow_teller(W, Z, W0, A)
        c = 1
        b = 5*A*c
        d = 0
        L = 0
        C = sf.shape_factor_gamow_teller(W, Z, W0, R, A, b, c, d, L)
    l = thecobs.Screening.screening_potential(Z)
    s = sf.atomic_screening(W, Z, R, l)

    sp = ph*f*l0*u*rc*rec*recCoul*C*s

    comb = np.stack((E, W, sp, ph, f, l0, rc, C, s, u, rec, recCoul), axis=1)

    df = pd.DataFrame(comb, columns = ['Energy', 'W', 'Spectrum', 'PhaseSpace', 'FermiFunction', 'L0', 'RadiativeCorrections', 'ShapeFactor', 'Screening', 'U', 'Recoil', 'CoulombRecoil'])

    return df

@st.cache
def convert_df(df):
    # IMPORTANT: Cache the conversion to prevent computation on every rerun
    return df.to_csv(index=False).encode('utf-8')

st.subheader('Transition info')

st.write("""
Current transition information:
""")

if sl_e0 > 0:
    zeff = z+1 if decay_type == 'Beta-' else -(z-1)
    df = calculateSpectrum(zeff, a, r, sl_e0, sl_e_step, beta_type)

    st.subheader('Electron spectrum')

    p = figure(title="Electron spectrum", x_axis_label='Kinetic energy [keV]', y_axis_label='dN/dE')

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
            pcorr = figure(title='Corrections', x_axis_label='Kinetic energy [keV]', y_axis_label='Correction', y_range=DataRange1d(only_visible=True))
            i = 0
            dfCorr = df.loc[:, ~df.columns.isin(['Energy', 'W', 'Spectrum', 'PhaseSpace'])]
            for column in dfCorr:
                pcorr.line(df['Energy'], df[column], legend_label=column, color=Category10[len(dfCorr.columns)][i], line_width=2)
                i+=1
            pcorr.legend.click_policy = "hide"

            st.bokeh_chart(pcorr, use_container_width=True)

        if st.checkbox('Show raw data'):
            df
else:
    st.error("Endpoint energy is less than 0. Can't calculate spectrum.")
