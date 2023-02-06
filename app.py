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

@st.cache
def downloadNubase():
    names = ['Z', 'A', 'halflife', 'halflifeUnc', 'spinParity']
    unitDict = {'My': 1e6*365.25*24*3600, 'ky': 1e3*365.25*24*3600, 'y': 365.25*24*3600, 'd': 24*3600., 'm': 60., 's': 1., 'ms': 1e-3, 'us': 1e-6, 'ns': 1e-9, 'ps': 1e-12, 'fs': 1e-15, 'as': 1e-18, 'zs': 1e-21, 'ys': 1e-24}
    data = []
    url = 'https://www-nds.iaea.org/amdc/ame2020/nubase_4.mas20.txt'
    r = requests.get(url).content
    f = StringIO(r.decode('utf-8'))
    skipLines = 25
    currentLine = 0
    for line in f:
        currentLine += 1
        if currentLine <= skipLines:
            continue
        else:
            a = int(line[0:3])
            zi = line[4:8]
            if zi[-1] != '0':
                continue
            z = int(zi[:-1])
            halflife = line[69:78]
            halflifeUnit = line[78:80].strip()
            halflifeUnc = line[81:88]
            spinParity = line[88:102].split(' ')[0].replace('#', '').replace('*', '')
            try:
                halflife = float(halflife)*unitDict[halflifeUnit]
                halflifeUnc = float(halflifeUnc)*unitDict[halflifeUnit]
            except:
                halflife = -1.
                halflifeUnc = 0.

            data.append([z, a, halflife, halflifeUnc, spinParity])
    return pd.DataFrame(data, columns=names)

apptitle = "The cob(s)web"

st.set_page_config(page_title=apptitle, page_icon=":eyeglasses:")

st.title("The Cob(s)web")

st.markdown("""
 Calculate beta spectrum shapes with state of the art methods
""")

st.sidebar.markdown("## Set the decay parameters")

st.sidebar.markdown("""Pulling data from Atomic Mass Evaluation 2020 and Nubase 2020 [link](https://www-nds.iaea.org/amdc/)
        """)

str_iso = st.sidebar.text_input('Isotope', value='1N')
z = 0.
a = 1.

m = re.search(r'\d+', str_iso)
if m:
    a = int(m.group(0))
    z = atoms.index(str(str_iso).replace(m.group(0), '').strip())
else:
    st.error('Not a valid isotope name. Expecting something like 6He or 45Ca.')

decay_type = st.sidebar.selectbox('Type of beta transition',
                                    ['Beta-', 'Beta+'])

halflife = 1.
Qvalue = 1000.
dfAME = downloadAME()
dfNubase = downloadNubase()

if decay_type == 'Beta-':
    Qvalue = dfAME.loc[(dfAME['Z'] == z) & (dfAME['A'] == a), 'mass'].values[0]-dfAME.loc[(dfAME['Z'] == (z+1)) & (dfAME['A'] == a), 'mass'].values[0]
else:
    Qvalue = dfAME.loc[(dfAME['Z'] == z) & (dfAME['A'] == a), 'mass'].values[0]-dfAME.loc[(dfAME['Z'] == (z-1)) & (dfAME['A'] == a), 'mass'].values[0]-2*ELECTRON_MASS_KEV

halflife = dfNubase.loc[(dfNubase['Z'] == z) & (dfNubase['A'] == a), 'halflife'].values[0]

with st.sidebar.expander("Nuclear data"):
    z = st.number_input('Proton number', min_value=0, max_value=120, step=1, value=z, help="Set the proton number of the initial state")
    a = st.number_input('Mass number', min_value=z, max_value=120, step=1, value=a, help="Set the mass number of the initial state")
    sl_r = st.number_input('Radius', min_value=0.01, max_value=100., step=0.1, value=approxNuclRadius(z, a), help="Set the rms nuclear radius in fm")
    sl_halflife = st.number_input('Partial falflife ($$t_{1/2}^{\\beta}$$)', min_value=0., max_value=1e32, value=halflife, step=1e-3, format='%e', help='Enter the partial halflife of the beta decay in seconds to calculate a log ft value')

    r = float(sl_r)*(5/3)**0.5*1e-15/NATURAL_LENGTH

    #str_iso.value = "%d%s" % (a, atoms[z])

sl_e0 = st.sidebar.number_input('Endpoint energy', min_value=0., max_value=20e3, value=Qvalue, step=1., help="Set the endpoint energy in keV")
sl_e_step = st.sidebar.number_input('Energy step', min_value=0.1, max_value=1000., value=1., step=1., help="Set the step energy in keV")

beta_type = st.sidebar.selectbox('Type of allowed transition',
                                    ['Gamow-Teller', 'Fermi', 'Mixed'])

mixing_ratio = 0.
if beta_type == 'Mixed':
    mixing_ratio_sl = st.sidebar.slider('Mixing ratio', -3., 3., -2.22, help='Ratio of Gamow-Teller and Fermi matrix elements for a mixed decay')
    mixing_ratio = float(mixing_ratio_sl)

bAc = 0.
dAc = 0.
Lambda = 0.
if beta_type == 'Gamow-Teller' or beta_type == 'Mixed':
    with st.sidebar.expander("Recoil-order matrix elements"):
        st.write("Specify the values for the recoil-order matrix elements following the notation in Hayen et al. Reviews of Modern Physics 90 (2018) 015008")
        bAc = st.number_input('b/Ac', min_value=-100., max_value=100.,value=5.0, step=0.1, help='Weak magnetism according to Holstein. Typical values are around 5.')
        dAc = st.number_input('d/Ac', min_value=-100., max_value=100., value=2., step=0.1, help='Induced tensor according to Holstein. Typical values in {-5, 5}.')
        Lambda = st.number_input(r'$$\Lambda$$', min_value=-100., max_value=100., value=5., step=0.1, help='Induced pseudoscalar according to Hayen et al. Typical values in {-5, 5}.')

@st.cache
def calculateSpectrum(Z, A, R, E0, E_step, beta_type, mixing_ratio=0, bAc=0, dAc=0, Lambda=0):
    E = np.arange(1., E0, E_step)

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
        b = bAc*A*c
        d = dAc*A*c
        L = Lambda
        C = sf.shape_factor_gamow_teller(W, Z, W0, R, A, b, c, d, L)
    elif beta_type == 'Mixed':
        norm = 1+mixing_ratio**2
        recF = sf.recoil_fermi(W, W0, A)
        recCoulF = sf.recoil_Coulomb_fermi(W, Z, W0, A)
        CF = sf.shape_factor_fermi(W, Z, W0, R)
        recGT = sf.recoil_gamow_teller(W, W0, A)
        recCoulGT = sf.recoil_Coulomb_gamow_teller(W, Z, W0, A)
        c = 1
        b = 5*A*c
        d = 0
        L = 0
        CGT = sf.shape_factor_gamow_teller(W, Z, W0, R, A, b, c, d, L)
        rec = 1+(recF-1+mixing_ratio**2*(recGT-1))/norm
        recCoul = 1+(recCoulF-1+mixing_ratio**2*(recCoulGT-1))/norm
        C = 1+(CF-1 + mixing_ratio**2*(CGT-1))/norm
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

if sl_e0 > 0:
    zeff = z+1 if decay_type == 'Beta-' else -(z-1)
    df = calculateSpectrum(zeff, a, r, sl_e0, sl_e_step, beta_type, mixing_ratio, bAc, dAc, Lambda)

    ftValuePS = np.sum(df['PhaseSpace'])*sl_e_step/ELECTRON_MASS_KEV*halflife
    ftValueFull = np.sum(df['Spectrum'])*sl_e_step/ELECTRON_MASS_KEV*halflife

    initState_str = "$$^{%d}$$%s" % (a, atoms[z])
    finalState_str = "$$^{%d}$$%s" % (a, atoms[abs(zeff)])

    st.subheader('Transition info')

    st.write("""
    Current transition information: %s   &rarr;   %s [$$E_0$$ = %.2f keV]

    Log ft values: 
    * Only phase space: %.5f
    * Full calculation: %.5f
    """ % (initState_str, finalState_str, sl_e0, np.log10(ftValuePS), np.log10(ftValueFull)))

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
    st.error("Endpoint energy is less than 0. Can't calculate spectrum. Did you mean to change the decay type(beta+/-)?")
