#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
calc_qi2009_udl.py

Calculate alpha/cluster-decay half-lives using Qi et al. (2009)
Universal Decay Law (UDL):

  log10(T1/2/s) = a*chi' + b*rho' + c
  chi' = Zc*Zd*sqrt(Ared/Qc)
  rho' = sqrt(Ared*Zc*Zd*(Ad^(1/3)+Ac^(1/3)))
  Ared = Ac*Ad/(Ac+Ad)

Default coefficient set: alpha fit from Qi et al. Table I:
  a=0.4065, b=-0.4311, c=-20.7889
Use --set alpha_fixed_a, cluster, cluster_fixed_a, all, all_fixed_a if desired.

Input CSV examples:
  Name,A,Z,Ealpha_MeV,Texp_s,Br_alpha,l
  273Ds_a,273,110,10.929,0.030,1.0,0

Optional cluster columns: Ac,Zc,Qc_MeV.
"""
from __future__ import annotations
import argparse, math
from typing import Optional
import pandas as pd

ALIASES={
 'A':['A','Ap','A_parent','Mass','MassNumber'], 'Z':['Z','Zp','Z_parent','Proton','ProtonNumber'],
 'N':['N','Np','N_parent','Neutron','NeutronNumber'], 'Ac':['Ac','A_cluster','A_c','Aemit','A_emitted'],
 'Zc':['Zc','Z_cluster','Z_c','Zemit','Z_emitted'],
 'Q':['Qc_MeV','Qcluster_MeV','Qalpha_MeV','Q_alpha_MeV','Qalpha','Q_alpha','Q_MeV','Q'],
 'E':['Ealpha_MeV','E_alpha_MeV','Ealpha','E_alpha','EXP','EXP_MeV','Ea','Ea_MeV'],
 'Texp':['Texp_s','T_exp_s','T0.5','T12_s','T1/2_s','half_life_s','tau_s'],
 'Br':['Br_alpha','Br','BR','Branch','b_alpha','branch_alpha'],
}
COEFFS={
 'alpha':        (0.4065, -0.4311, -20.7889),
 'cluster':      (0.3671, -0.3296, -26.2681),
 'all':          (0.3949, -0.3693, -23.7615),
 'alpha_fixed_a':(0.4314, -0.4608, -21.9453),
 'cluster_fixed_a':(0.4314, -0.3921, -32.7044),
 'all_fixed_a':  (0.4314, -0.4087, -25.7725),
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
def qi_udl_log10T(A:int,Z:int,Qc:float,Ac:int=4,Zc:int=2,coeff_set='alpha'):
    if Qc<=0: raise ValueError('Qc must be positive')
    Ad=A-Ac; Zd=Z-Zc
    if Ad<=0 or Zd<=0: raise ValueError('Invalid daughter nucleus')
    Ared=Ac*Ad/(Ac+Ad)
    chip=Zc*Zd*math.sqrt(Ared/Qc)
    rhop=math.sqrt(Ared*Zc*Zd*(Ad**(1/3)+Ac**(1/3)))
    a,b,c=COEFFS[coeff_set]
    return a*chip+b*rhop+c, chip, rhop

def main():
    ap=argparse.ArgumentParser(description='Qi et al. 2009 UDL half-life calculator')
    ap.add_argument('input_csv'); ap.add_argument('-o','--output',default='output_qi2009_udl.csv')
    ap.add_argument('--set',choices=list(COEFFS),default='alpha',help='Coefficient set from Qi Table I')
    ap.add_argument('--exp-is-qalpha',action='store_true',help='Treat Ealpha/EXP as Q directly')
    args=ap.parse_args()
    df=pd.read_csv(args.input_csv); df.columns=[str(c).strip().replace('\ufeff','') for c in df.columns]
    cA,cZ,cAc,cZc,cQ,cE,cT,cBr=[find_col(df,k) for k in ['A','Z','Ac','Zc','Q','E','Texp','Br']]
    if cA is None or cZ is None: raise ValueError(f'Need A and Z. Found {list(df.columns)}')
    if cQ is None and cE is None: raise ValueError(f'Need Qalpha/Qc or Ealpha. Found {list(df.columns)}')
    rows=[]
    for i,row in df.iterrows():
        A=int(round(fnum(row[cA]))); Z=int(round(fnum(row[cZ])))
        Ac=int(round(fnum(row[cAc],4))) if cAc else 4; Zc=int(round(fnum(row[cZc],2))) if cZc else 2
        if cQ: Q=fnum(row[cQ]); qsrc=cQ
        else:
            E=fnum(row[cE]);
            if 'kev' in norm(cE): E/=1000.0
            Q=E if args.exp_is_qalpha else E*A/(A-Ac); qsrc=(cE+' treated as Q' if args.exp_is_qalpha else cE+f' converted by A/(A-{Ac})')
        logT,chip,rhop=qi_udl_log10T(A,Z,Q,Ac,Zc,args.set)
        out={'Ad':A-Ac,'Zd':Z-Zc,'Ac_used':Ac,'Zc_used':Zc,'Qc_QiUDL_MeV':Q,'Q_source':qsrc,
             'coeff_set_QiUDL':args.set,'chi_prime':chip,'rho_prime':rhop,'log10T_Qi2009_UDL_s':logT,'T_Qi2009_UDL_s':10**logT}
        if cT:
            T=fnum(row[cT]); Br=fnum(row[cBr],1.0) if cBr else 1.0
            out['Talpha_exp_s']=T/Br; out['HF_Qi2009_UDL']=(T/Br)/(10**logT); out['log10HF_Qi2009_UDL']=math.log10(out['HF_Qi2009_UDL'])
        rows.append(out)
    pd.concat([df.reset_index(drop=True),pd.DataFrame(rows)],axis=1).to_csv(args.output,index=False)
    print(f'Wrote {args.output}')
if __name__=='__main__': main()
