#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
calc_xu2022_unified.py

Calculate alpha-decay half-lives using Xu et al. (2022) unified formula:

 log10(T1/2/s) = F(Z)*(Ad/(Ap*Qalpha))^1/2 * [acos(sqrt(X))-sqrt(X*(1-X))]
                 - 20.446 + C(Z,N) + d*l(l+1) + h

 F(Z) = 28.274*sqrt(Z) + 2920.347/Z - 204.086
 X = r0*(Ad^(1/3)+Aalpha^(1/3))*Qalpha/(2*Zd*e2)
 r0 = 1.2249 fm, e2 = 1.4399764 MeV fm, d = 0.0669
 h = 0 for even-even, 0.2018 for odd-A, 0.4036 for odd-odd parent nuclei.

C(Z,N) is defined by Xu et al. only for two regions near Pb:
  78<=Z<=82 and 100<=N<126, or 82<Z<=90 and 110<=N<=126.
Outside those ranges this script uses C=0 and marks C_status='outside_defined_region'.

Input CSV examples:
  Name,A,Z,Ealpha_MeV,Texp_s,Br_alpha,l
  273Ds_a,273,110,10.929,0.030,1.0,0
"""
from __future__ import annotations
import argparse, math
from typing import Optional
import pandas as pd

R0=1.2249
E2=1.4399764
D_L=0.0669
AALPHA=4
ZALPHA=2
ALIASES={
 'A':['A','Ap','A_parent','Mass','MassNumber'], 'Z':['Z','Zp','Z_parent','Proton','ProtonNumber'],
 'N':['N','Np','N_parent','Neutron','NeutronNumber'],
 'Q':['Qalpha_MeV','Q_alpha_MeV','Qalpha','Q_alpha','Q_MeV','Q'],
 'E':['Ealpha_MeV','E_alpha_MeV','Ealpha','E_alpha','EXP','EXP_MeV','Ea','Ea_MeV'],
 'l':['l','L','ell','DeltaL','Delta_L','lmin','l_min'],
 'Texp':['Texp_s','T_exp_s','T0.5','T12_s','T1/2_s','half_life_s','tau_s'],
 'Br':['Br_alpha','Br','BR','Branch','b_alpha','branch_alpha'],
}
def norm(s): return str(s).strip().replace('\ufeff','').lower().replace(' ','').replace('-','_')
def find_col(df,key)->Optional[str]:
    mp={norm(c):c for c in df.columns}
    for a in ALIASES[key]:
        if norm(a) in mp: return mp[norm(a)]
    return None
def fnum(x, default=math.nan):
    if pd.isna(x): return default
    if isinstance(x,str) and x.strip()=='': return default
    return float(x)
def F_of_Z(Z): return 28.274*math.sqrt(Z)+2920.347/Z-204.086
def C_shell(Z,N):
    if 78 <= Z <= 82 and 100 <= N < 126:
        return 1.547 - 0.077*(82-Z) - 0.050*(126-N), 'defined_region_1'
    if 82 < Z <= 90 and 110 <= N <= 126:
        return 1.397 - 0.116*(Z-82) - 0.061*(126-N), 'defined_region_2'
    return 0.0, 'outside_defined_region'
def h_blocking(Z,N):
    if Z%2==0 and N%2==0: return 0.0, 'even-even'
    if Z%2==1 and N%2==1: return 0.4036, 'odd-odd'
    return 0.2018, 'odd-A'
def xu_log10T(A:int,Z:int,Q:float,ell:float=0.0):
    if Q<=0: raise ValueError('Qalpha must be positive')
    Ad=A-AALPHA; Zd=Z-ZALPHA; N=A-Z
    X=R0*(Ad**(1/3)+AALPHA**(1/3))*Q/(2*Zd*E2)
    if not (0 < X < 1): raise ValueError(f'Invalid X={X}; need 0<X<1')
    bracket=math.acos(math.sqrt(X))-math.sqrt(X*(1-X))
    C,cs=C_shell(Z,N); h,case=h_blocking(Z,N)
    logT=F_of_Z(Z)*math.sqrt(Ad/(A*Q))*bracket - 20.446 + C + D_L*ell*(ell+1.0) + h
    return logT, X, bracket, C, cs, h, case

def main():
    ap=argparse.ArgumentParser(description='Xu et al. 2022 unified alpha half-life calculator')
    ap.add_argument('input_csv'); ap.add_argument('-o','--output',default='output_xu2022_unified.csv')
    ap.add_argument('--exp-is-qalpha',action='store_true',help='Treat Ealpha/EXP column as Qalpha directly')
    args=ap.parse_args()
    df=pd.read_csv(args.input_csv); df.columns=[str(c).strip().replace('\ufeff','') for c in df.columns]
    cA,cZ,cN,cQ,cE,cl,cT,cBr=[find_col(df,k) for k in ['A','Z','N','Q','E','l','Texp','Br']]
    if cA is None or cZ is None: raise ValueError(f'Need A and Z columns. Found {list(df.columns)}')
    if cQ is None and cE is None: raise ValueError(f'Need Qalpha or Ealpha column. Found {list(df.columns)}')
    rows=[]
    for i,row in df.iterrows():
        A=int(round(fnum(row[cA]))); Z=int(round(fnum(row[cZ]))); N=int(round(fnum(row[cN]))) if cN else A-Z
        if cQ: Q=fnum(row[cQ]); qsrc=cQ
        else:
            E=fnum(row[cE]);
            if 'kev' in norm(cE): E/=1000.0
            Q=E if args.exp_is_qalpha else E*A/(A-4.0); qsrc=(cE+' treated as Qalpha' if args.exp_is_qalpha else cE+' converted by A/(A-4)')
        ell=fnum(row[cl],0.0) if cl else 0.0
        logT,X,br,C,cs,h,case=xu_log10T(A,Z,Q,ell)
        out={'N_calc':N,'Ad':A-4,'Zd':Z-2,'Qalpha_Xu_MeV':Q,'Q_source':qsrc,'l_used':ell,
             'X_Xu':X,'bracket_Xu':br,'F_Z_Xu':F_of_Z(Z),'C_ZN_Xu':C,'C_status':cs,'h_blocking_Xu':h,'parity_case':case,
             'log10T_Xu2022_s':logT,'T_Xu2022_s':10**logT}
        if cT:
            T=fnum(row[cT]); Br_alpha=fnum(row[cBr],1.0) if cBr else 1.0
            out['Talpha_exp_s']=T/Br_alpha; out['HF_Xu2022']=(T/Br_alpha)/(10**logT); out['log10HF_Xu2022']=math.log10(out['HF_Xu2022'])
        rows.append(out)
    pd.concat([df.reset_index(drop=True),pd.DataFrame(rows)],axis=1).to_csv(args.output,index=False)
    print(f'Wrote {args.output}')
if __name__=='__main__': main()
