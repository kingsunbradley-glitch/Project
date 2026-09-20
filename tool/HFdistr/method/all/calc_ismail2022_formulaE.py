#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
calc_ismail2022_formulaE.py

Calculate alpha/cluster-decay half-lives using Ismail et al. (2022)
Formula E (improved NRDX with angular momentum, isospin and parity):

  log10(T1/2/s) = a*sqrt(mu)*Zc*Zd/sqrt(Qc)
                 + b*sqrt(mu)*sqrt(Zc*Zd) + c
                 + d*l(l+1) + e*I + f*I^2 + g*(1 - (-1)^l)

where mu = Ac*Ad/(Ac+Ad), I=(N_parent-Z_parent)/A_parent.
For ordinary alpha decay, Ac=4 and Zc=2 are used by default.

Input CSV examples:
  Name,A,Z,Ealpha_MeV,Texp_s,Br_alpha,l
  273Ds_a,273,110,10.929,0.030,1.0,0

Optional columns for cluster decay:
  Ac,Zc,Qc_MeV

Output: input columns plus Qc, log10T, Tcalc and HF if Texp_s exists.
"""

from __future__ import annotations
import argparse, math
from pathlib import Path
from typing import Optional
import pandas as pd

ALIASES = {
    'A': ['A','Ap','A_parent','Mass','MassNumber'],
    'Z': ['Z','Zp','Z_parent','Proton','ProtonNumber'],
    'N': ['N','Np','N_parent','Neutron','NeutronNumber'],
    'Ac': ['Ac','A_cluster','A_c','Aemit','A_emitted'],
    'Zc': ['Zc','Z_cluster','Z_c','Zemit','Z_emitted'],
    'Q': ['Qc_MeV','Qcluster_MeV','Qalpha_MeV','Q_alpha_MeV','Qalpha','Q_alpha','Q_MeV','Q'],
    'E': ['Ealpha_MeV','E_alpha_MeV','Ealpha','E_alpha','EXP','EXP_MeV','Ea','Ea_MeV'],
    'l': ['l','L','ell','DeltaL','Delta_L','lmin','l_min'],
    'Texp': ['Texp_s','T_exp_s','T0.5','T12_s','T1/2_s','half_life_s','tau_s'],
    'Br': ['Br_alpha','Br','BR','Branch','b_alpha','branch_alpha'],
}

# Table 5 in Ismail et al. 2022, Formula E. Keys classify parent Z,N parity.
PARAMS = {
    'even-even': dict(a=0.410918, b=-1.42120, c=-15.2909, d=0.0,      e=10.73868, f=-56.9910, g=0.0),
    'even-odd':  dict(a=0.410796, b=-1.42969, c=-14.2146, d=0.022595, e=-2.25920, f=-0.53224, g=0.502634),
    'odd-even':  dict(a=0.424757, b=-1.31517, c=-17.6008, d=0.020609, e=-12.09466,f=10.0859,  g=-0.134478),
    'odd-odd':   dict(a=0.421951, b=-1.39998, c=-16.3099, d=0.030543, e=0.96066,  f=-12.2475, g=0.227788),
    'all':       dict(a=0.408505, b=-1.41695, c=-14.6766, d=0.027154, e=4.13841,  f=-24.9919, g=0.244732),
}

def norm(s: str) -> str:
    return str(s).strip().replace('\ufeff','').lower().replace(' ','').replace('-','_')

def find_col(df: pd.DataFrame, key: str) -> Optional[str]:
    mp = {norm(c): c for c in df.columns}
    for a in ALIASES[key]:
        if norm(a) in mp:
            return mp[norm(a)]
    return None

def fnum(x, default=math.nan) -> float:
    if pd.isna(x): return default
    if isinstance(x,str) and x.strip()=='': return default
    return float(x)

def parity_case(Z:int, N:int) -> str:
    if Z%2==0 and N%2==0: return 'even-even'
    if Z%2==0 and N%2==1: return 'even-odd'
    if Z%2==1 and N%2==0: return 'odd-even'
    return 'odd-odd'

def ismail_log10T(A:int, Z:int, Qc:float, ell:float=0.0, Ac:int=4, Zc:int=2, param_set: str|None=None) -> float:
    if Qc <= 0: raise ValueError('Qc must be positive')
    Ad, Zd = A - Ac, Z - Zc
    if Ad <= 0 or Zd <= 0: raise ValueError('Invalid daughter nucleus from A,Z,Ac,Zc')
    N = A - Z
    key = param_set or parity_case(Z, N)
    p = PARAMS[key]
    mu = Ac * Ad / (Ac + Ad)
    I = (N - Z) / A
    return (p['a']*math.sqrt(mu)*Zc*Zd/math.sqrt(Qc)
            + p['b']*math.sqrt(mu)*math.sqrt(Zc*Zd) + p['c']
            + p['d']*ell*(ell+1.0) + p['e']*I + p['f']*I*I
            + p['g']*(1.0 - ((-1.0)**int(round(ell)))))

def main():
    ap = argparse.ArgumentParser(description='Ismail et al. 2022 Formula E half-life calculator')
    ap.add_argument('input_csv')
    ap.add_argument('-o','--output',default='output_ismail2022_formulaE.csv')
    ap.add_argument('--exp-is-qalpha',action='store_true',help='Treat Ealpha/EXP column as Qalpha/Qc directly')
    ap.add_argument('--param-set',choices=list(PARAMS),default=None,help='Override parity-specific coefficients')
    args = ap.parse_args()
    df = pd.read_csv(args.input_csv)
    df.columns = [str(c).strip().replace('\ufeff','') for c in df.columns]
    cA,cZ,cN,cAc,cZc,cQ,cE,cl,cT,cBr = [find_col(df,k) for k in ['A','Z','N','Ac','Zc','Q','E','l','Texp','Br']]
    if cA is None or cZ is None: raise ValueError(f'Need A and Z columns. Found {list(df.columns)}')
    if cQ is None and cE is None: raise ValueError(f'Need Qalpha/Qc or Ealpha column. Found {list(df.columns)}')
    rows=[]
    for i,row in df.iterrows():
        A=int(round(fnum(row[cA]))); Z=int(round(fnum(row[cZ]))); N=int(round(fnum(row[cN]))) if cN else A-Z
        Ac=int(round(fnum(row[cAc],4))) if cAc else 4
        Zc=int(round(fnum(row[cZc],2))) if cZc else 2
        if cQ:
            Q=fnum(row[cQ]); qsrc=cQ
        else:
            E=fnum(row[cE]);
            if 'kev' in norm(cE): E/=1000.0
            Q=E if args.exp_is_qalpha else E*A/(A-Ac)
            qsrc=(cE+' treated as Q' if args.exp_is_qalpha else cE+f' converted by A/(A-{Ac})')
        ell=fnum(row[cl],0.0) if cl else 0.0
        case=parity_case(Z,N)
        logT=ismail_log10T(A,Z,Q,ell,Ac,Zc,args.param_set)
        out={'N_calc':N,'Ad':A-Ac,'Zd':Z-Zc,'Ac_used':Ac,'Zc_used':Zc,'Qc_Ismail_MeV':Q,
             'Q_source':qsrc,'parity_case':case,'param_set_used':args.param_set or case,
             'l_used':ell,'log10T_Ismail2022_FE_s':logT,'T_Ismail2022_FE_s':10**logT}
        if cT:
            T=fnum(row[cT]); Br=fnum(row[cBr],1.0) if cBr else 1.0
            out['Talpha_exp_s']=T/Br; out['HF_Ismail2022_FE']=(T/Br)/(10**logT); out['log10HF_Ismail2022_FE']=math.log10(out['HF_Ismail2022_FE'])
        rows.append(out)
    pd.concat([df.reset_index(drop=True),pd.DataFrame(rows)],axis=1).to_csv(args.output,index=False)
    print(f'Wrote {args.output}')

if __name__ == '__main__': main()
