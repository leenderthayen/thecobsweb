import streamlit as st
import csv
import numpy as np
import pandas as pd
import re
import thecobs.SpectralFunctions as sf
from thecobs.Constants import *
import thecobs.Functions as fu
import thecobs.Screening
import thecobs.Generator as g

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
    names = ['Z', 'A', 'halflife', 'halflifeUnc', 'spinParity', 'branchingB-', 'branchingB+']
    unitDict = {'Qy': 1e30*365.25*24*3600, 'Ry': 1e27*365.25*24*3600, 'Yy': 1e24*365.25*24*3600, 'Zy': 1e21*365.25*24*3600, 'Ey': 1e18*365.25*24*3600,'Py': 1e15*365.25*24*3600, 'Ty': 1e12*365.25*24*3600, 'Gy': 1e9*365.25*24*3600, 'My': 1e6*365.25*24*3600, 'ky': 1e3*365.25*24*3600, 'y': 365.25*24*3600, 'd': 24*3600., 'h': 3600., 'm': 60., 's': 1., 'ms': 1e-3, 'us': 1e-6, 'ns': 1e-9, 'ps': 1e-12, 'fs': 1e-15, 'as': 1e-18, 'zs': 1e-21, 'ys': 1e-24}
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
            branches = line[119:209]
            branchingRatioBm = 100.
            branchingRatioBp = 100.
            try:
                halflife = float(halflife)*unitDict[halflifeUnit]
                halflifeUnc = float(halflifeUnc)*unitDict[halflifeUnit]
            except:
                halflife = -1.
                halflifeUnc = 0.
            try:
                m = re.search(r'B-=\d+[\.\d+]*', branches)
                if m:
                    branchingRatioBm = float(m.group(0).split('B-=')[1])
                m = re.search(r'B+=\d+[\.\d+]*', branches)
                if m:
                    branchingRatioBp = float(m.group(0).split('B+=')[1])
            except:
                pass

            data.append([z, a, halflife, halflifeUnc, spinParity, branchingRatioBm, branchingRatioBp])
    return pd.DataFrame(data, columns=names)

apptitle = "The cob(s)web"

st.set_page_config(page_title=apptitle, page_icon=":eyeglasses:")

st.title("The Cob(s)web")

st.markdown("""
The calculation of beta spectra, now on the web!

State of the art methods as described [here](#citation-and-about).
""")

st.sidebar.markdown("## Decay parameters")

str_iso = st.sidebar.text_input('Isotope', value='1N', help="Enter the initial state as an isotope name, like 19Ne or 63Ni")
z = 0.
a = 1.

m = re.search(r'\d+', str_iso)
if m:
    a = int(m.group(0))
    z = atoms.index(str(str_iso).replace(m.group(0), '').strip().capitalize())
else:
    st.error('Not a valid isotope name. Expecting something like 6He or 45Ca.')

decay_type = st.sidebar.selectbox('Type of beta transition',
                                    ['Beta-', 'Beta+'])

halflife = 1.
branchingRatio = 100.
Qvalue = 1000.
spinParityInit = '1/2+'
spinParityFinal = '1/2+'
dfAME = downloadAME()
dfNubase = downloadNubase()

halflife = dfNubase.loc[(dfNubase['Z'] == z) & (dfNubase['A'] == a), 'halflife'].values[0]
spinParityInit = dfNubase.loc[(dfNubase['Z'] == z) & (dfNubase['A'] == a), 'spinParity'].values[0]

if decay_type == 'Beta-':
    Qvalue = dfAME.loc[(dfAME['Z'] == z) & (dfAME['A'] == a), 'mass'].values[0]-dfAME.loc[(dfAME['Z'] == (z+1)) & (dfAME['A'] == a), 'mass'].values[0]
    spinParityFinal = dfNubase.loc[(dfNubase['Z'] == z+1) & (dfNubase['A'] == a), 'spinParity'].values[0]
    branchingRatio = dfNubase.loc[(dfNubase['Z'] == z) & (dfNubase['A'] == a), 'branchingB-'].values[0]
else:
    Qvalue = dfAME.loc[(dfAME['Z'] == z) & (dfAME['A'] == a), 'mass'].values[0]-dfAME.loc[(dfAME['Z'] == (z-1)) & (dfAME['A'] == a), 'mass'].values[0]-2*ELECTRON_MASS_KEV
    spinParityFinal = dfNubase.loc[(dfNubase['Z'] == z-1) & (dfNubase['A'] == a), 'spinParity'].values[0]
    branchingRatio = dfNubase.loc[(dfNubase['Z'] == z) & (dfNubase['A'] == a), 'branchingB+'].values[0]

with st.sidebar.expander("Nuclear data"):
    z = st.number_input('Proton number', min_value=0, max_value=120, step=1, value=z, help="Set the proton number of the initial state")
    a = st.number_input('Mass number', min_value=z, max_value=300, step=1, value=a, help="Set the mass number of the initial state")
    sl_r = st.number_input('Radius', min_value=0.005, max_value=100., step=0.1, value=approxNuclRadius(z, a), help="Set the rms nuclear radius in fm")
    sl_halflife = st.number_input('Halflife ($$t_{1/2}^{\\beta}$$)', min_value=0., max_value=1e32, value=halflife, step=1e-3, format='%e', help='Enter the halflife of the initial state in seconds to calculate a log ft value')
    sl_branching = st.number_input('Branching ratio (%)', min_value=1e-3, max_value=100., value=branchingRatio, step=1., help="Enter the beta decay branching ratio as a percentage")
    str_spinParityInit = st.text_input('Initial spin-parity', value=spinParityInit, help="Define the spin-parity of the initial state, like 1/2+, 7/2-, or 2+")
    str_spinParityFinal = st.text_input('Final spin-parity', value=spinParityFinal, help="Define the spin-parity of the final state, like 1/2+, 7/2-, or 2+")

    r = float(sl_r)*(5/3)**0.5*1e-15/NATURAL_LENGTH
    spinParityInit = str(str_spinParityInit)
    spinParityFinal = str(str_spinParityFinal)
    branchingRatio = float(sl_branching)

    #str_iso.value = "%d%s" % (a, atoms[z])

sl_e0 = st.sidebar.number_input('Endpoint energy', min_value=0., max_value=20e3, value=Qvalue, step=1., help="Set the endpoint energy in keV")
sl_e_step = st.sidebar.number_input('Energy step', min_value=0.1, max_value=1000., value=0.1, step=1., help="Set the step energy in keV")

dJ, f, u = fu.determineForbiddenness(spinParityInit, spinParityFinal)
transIndex = 0
if f != 0:
  transIndex += f+2

beta_type = st.sidebar.selectbox('Type of transition',
        ['A: Gamow-Teller', 'A: Fermi', 'A: Mixed', '1FU', '2FU', '3FU', '4FU'], index=transIndex, help="Choose between the different types of allowed (A) or uniquely n-th forbidden (nFU) transitions.")

if 'FU' in beta_type:
    f = int(beta_type[0])
    dJ = f+1
    u = True

mixing_ratio = 0.
if beta_type == 'A: Mixed':
    mixing_ratio_sl = st.sidebar.slider('Mixing ratio', -3., 3., -2.22, help='Ratio of Gamow-Teller and Fermi matrix elements for a mixed decay')
    mixing_ratio = float(mixing_ratio_sl)

bAc = 0.
dAc = 0.
Lambda = 0.
if beta_type == 'A: Gamow-Teller' or beta_type == 'A: Mixed':
    with st.sidebar.expander("Recoil-order matrix elements"):
        st.write("Specify the values for the recoil-order matrix elements following the notation in Hayen et al. [Reviews of Modern Physics 90 (2018) 015008](https://journals.aps.org/rmp/abstract/10.1103/RevModPhys.90.015008)")
        bAc = st.number_input('b/Ac', min_value=-100., max_value=100.,value=5.0, step=0.1, help='Weak magnetism according to Holstein. Typical values are around 5.')
        dAc = st.number_input('d/Ac', min_value=-100., max_value=100., value=2., step=0.1, help='Induced tensor according to Holstein. Typical values in {-5, 5}.')
        Lambda = st.number_input(r'$$\Lambda$$', min_value=-100., max_value=100., value=5., step=0.1, help='Induced pseudoscalar according to Hayen et al. Typical values in {-5, 5}.')

st.sidebar.markdown("""Extracting data from [Atomic Mass Evaluation 2020 and Nubase 2020](https://www-nds.iaea.org/amdc/)
        """)
st.sidebar.info("Q values derived from AME2020 neglect electron binding energy differences and can differ from [ENSDF/NuDat](https://www.nndc.bnl.gov/nudat3/)", icon="ℹ️")

@st.cache
def calculateSpectrum(Z, A, R, E0, E_step, dJ, beta_type, mixing_ratio=0, bAc=0, dAc=0, Lambda=0):
    E = np.arange(1., E0, E_step)

    W0 = 1 + E0/ELECTRON_MASS_KEV

    W = 1 + E/ELECTRON_MASS_KEV

    ph = sf.phase_space(W, W0)
    f = sf.fermi_function(W, Z, R)
    l0 = sf.finite_size_L0(W, Z, R)
    u = sf.finite_size_U_fermi(W, Z)
    rc = sf.radiative_correction(W, Z, W0, R)
    rec = np.ones(len(E))
    recCoul = np.ones(len(E))
    C = np.ones(len(E))

    if beta_type == 'A: Fermi':
        rec = sf.recoil_fermi(W, W0, A)
        recCoul = sf.recoil_Coulomb_fermi(W, Z, W0, A)
        C = sf.shape_factor_fermi(W, Z, W0, R)
    elif beta_type == 'A: Gamow-Teller':
        rec = sf.recoil_gamow_teller(W, W0, A)
        recCoul = sf.recoil_Coulomb_gamow_teller(W, Z, W0, A)
        c = 1
        b = bAc*A*c
        d = dAc*A*c
        L = Lambda
        C = sf.shape_factor_gamow_teller(W, Z, W0, R, A, b, c, d, L)
    elif beta_type == 'A: Mixed':
        norm = 1+mixing_ratio**2
        recF = sf.recoil_fermi(W, W0, A)
        recCoulF = sf.recoil_Coulomb_fermi(W, Z, W0, A)
        CF = sf.shape_factor_fermi(W, Z, W0, R)
        recGT = sf.recoil_gamow_teller(W, W0, A)
        recCoulGT = sf.recoil_Coulomb_gamow_teller(W, Z, W0, A)
        c = 1
        b = bAc*A*c
        d = dAc*A*c
        L = Lambda
        CGT = sf.shape_factor_gamow_teller(W, Z, W0, R, A, b, c, d, L)
        rec = 1+(recF-1+mixing_ratio**2*(recGT-1))/norm
        recCoul = 1+(recCoulF-1+mixing_ratio**2*(recCoulGT-1))/norm
        C = 1+(CF-1 + mixing_ratio**2*(CGT-1))/norm
    elif 'FU' in beta_type:
        C = sf.shape_factor_unique_forbidden(W, dJ, W0, Z, R)
    l = thecobs.Screening.screening_potential(Z)
    s = sf.atomic_screening(W, Z, R, l)
    X = np.ones(len(W))
    if Z > 0:
        exParsSim = g.getExchangeParamsSimkovic(Z)
        X = sf.atomic_exchange_simkovic(W, exParsSim)

    sp = ph*f*l0*u*rc*rec*recCoul*C*s*X

    comb = np.stack((E, W, sp, ph, f, l0, rc, C, s, X, u, rec, recCoul), axis=1)

    df = pd.DataFrame(comb, columns = ['Energy', 'W', 'Spectrum', 'PhaseSpace', 'FermiFunction', 'L0', 'RadiativeCorrections', 'ShapeFactor', 'Screening', 'Exchange', 'U', 'Recoil', 'CoulombRecoil'])

    return df

@st.cache
def convert_df(df):
    # IMPORTANT: Cache the conversion to prevent computation on every rerun
    return df.to_csv(index=False).encode('utf-8')

if sl_e0 > 0:
    zeff = z+1 if decay_type == 'Beta-' else -(z-1)
    df = calculateSpectrum(zeff, a, r, sl_e0, sl_e_step, dJ, beta_type, mixing_ratio, bAc, dAc, Lambda)

    initState_str = "$$^{%d}$$%s [%s]" % (a, atoms[z], spinParityInit)
    finalState_str = "$$^{%d}$$%s [%s]" % (a, atoms[abs(zeff)], spinParityFinal)

    st.subheader('Transition info')

    st.write("""
    Current transition information: %s &nbsp;&nbsp; &rarr; &nbsp;&nbsp; %s &nbsp;&nbsp; \n\r$$E_0$$ = %.2f keV &nbsp;&nbsp; BR=%.2f %s &nbsp;&nbsp; $$t_{1/2}$$ = %s s
    """ % (initState_str, finalState_str, sl_e0, branchingRatio, '(%)', 'stable' if halflife == -1 else '%.2e' % halflife))

    effHalflife = halflife/(branchingRatio/100.)

    if halflife < 0:
        st.warning("No halflife found. Setting $$t_{1/2}$$ = 1 s for calculating log ft.")
        effHalflife = 1.

    ftValuePS = np.sum(df['PhaseSpace'])*sl_e_step/ELECTRON_MASS_KEV*effHalflife
    ftValueFPS = np.sum(df['PhaseSpace']*df['FermiFunction'])*sl_e_step/ELECTRON_MASS_KEV*effHalflife
    ftValueFull = np.sum(df['Spectrum'])*sl_e_step/ELECTRON_MASS_KEV*effHalflife

    st.write("""
    Log ft values: [:question:](#comparative-halflife)
    * Only phase space: %.5f
    * With Fermi function: %.5f
    * Full calculation: %.5f
    """ % (np.log10(ftValuePS), np.log10(ftValueFPS), np.log10(ftValueFull)))

    st.subheader('Electron spectrum')

    normPS = normFPS = normFull = 1.

    if st.checkbox('Show normalized spectra'):
        normPS = 1./(ftValuePS/effHalflife)
        normFPS = 1./(ftValueFPS/effHalflife)
        normFull = 1./(ftValueFull/effHalflife)

    p = figure(title='Beta spectrum', x_axis_label='Kinetic energy [keV]', y_axis_label='dN/dE', y_range=DataRange1d(only_visible=True))

    p.line(df['Energy'], df['Spectrum']*normFull, legend_label='All corrections', line_width=2, muted_alpha=0.2, color=Category10[4][0])
    p.line(df['Energy'], df['PhaseSpace']*df['FermiFunction']*normFPS, legend_label='Fermi function', line_width=2, muted_alpha=0.2, color=Category10[4][1])
    p.line(df['Energy'], df['PhaseSpace']*normPS, legend_label='Phase Space', line_width=2, muted_alpha=0.2, color=Category10[4][3])

    p.legend.click_policy = "hide"

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

st.subheader("Comparative halflife")

st.write("""
The comparative half-life, or $$ft$$ value, describes the overlap between initial and final states involved in the $$\\beta$$ decay and is defined as follows

$$
ft = \\frac{K}{G_F^2V_{ud}^2}\\frac{1}{|\langle f | H_{\\beta}(0)|i\\rangle|^2}
$$
with $$K = 8120.27648(26) \cdot 10^{-10}$$ GeV$$^{-4}$$s. As the matrix element can vary by several orders of magnitude, it is convenient and conventional to instead report the base 10 logarithm of the $$ft$$ value as was done above.""")

with st.expander("More information"):
    st.write("""It does so by removing the effects from phase space, better known as the integral of the $$\\beta$$ spectrum or $$f$$ value,

$$
f = \int_1^{W_0} ~pW(W-W_0)^2 F(Z, W) C(Z, W) \ldots dW
$$
where $$W = 1+E_{kin}/m_ec^2 (p = \sqrt{W^2-1})$$ are the $$\\beta$$ particle dimensionless total energy (momentum) using $$m_e=c=\hbar=1$$, and combining it with an experimentally measured half-life, $$t_{1/2}$$. 

Depending on the type of transition, typical values are distributed as shown in the plot below from [these lecture notes](https://link.springer.com/chapter/10.1007/978-3-540-85839-3_4)
        """)

    st.image("https://www.researchgate.net/profile/Berta-Rubio/publication/225687961/figure/fig4/AS:302550210367490@1449144999109/Numbers-of-known-allowed-upper-panel-and-forbidden-lower-panel-transitions-of.png")

st.subheader("Citation and about")

st.write("""
The calculation of beta spectra uses thecobs([Github](https://github.com/leenderthayen/thecobs) and [Pypi](https://pypi.org/project/thecobs/)), a python port of the BSG software written in C++ ([Github](https://github.com/leenderthayen/BSG) and [documentation](https://bsg.readthedocs.io/en/latest/)).

Written by Leendert Hayen. If you use these results in a scientific publication, please cite

* L. Hayen et al., [Reviews of Modern Physics 90 (2018) 015008](https://journals.aps.org/rmp/abstract/10.1103/RevModPhys.90.015008)
* L. Hayen and N. Severijns, [Computer Physics Communications 240 (2019) 152](https://www.sciencedirect.com/science/article/abs/pii/S0010465519300645)
* O. Nitescu et al., [Physical Review C 107 (2023) 025501](https://journals.aps.org/prc/pdf/10.1103/PhysRevC.107.025501)

Find me on [Google Scholar](https://scholar.google.com/citations?user=2TlxGBkAAAAJ&hl=en).
        """)
