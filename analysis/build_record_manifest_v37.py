from pathlib import Path
import pandas as pd
ROOT=Path(r'E:\KimiData\kimi\Workspaces\论文写作')
old=ROOT/'work_v36'/'nhanes_rerun_v36_complete'
out=ROOT/'work_topic2_v37'/'results'
reg=pd.read_csv(old/'nhanes_registry_table.csv')
# GLM7-b is a second record for the same expression, not a 26th expression.
base=reg.loc[reg['index'].eq('GLM7')].iloc[0].copy()
base['index']='GLM7-b'; base['analysis_role']='dual_base_sensitivity'
base['component_terms']='sex;log10(age);log10(BMI);metab'
base['reference_class']='rank_preserving_exact'
reg=pd.concat([reg,pd.DataFrame([base])],ignore_index=True)
flow=pd.read_csv(old/'nhanes_record_flow.csv')
flow=flow.rename(columns={'record':'index'})
reg=reg.merge(flow[['index','role','evaluated_in']],on='index',how='left',validate='one_to_one')
rows=[]
for r in reg.itertuples(index=False):
    main=(r.analysis_role=='primary_candidate')
    arm_a=main
    # GLM7-b is the full embedded-base sensitivity record, not an Arm A main record.
    arm_b=r.index in ['GLM7','TyG','TyG_BMI','TyG_WC','TyG_WHtR','METS_IR','LAP','VAI','CVAI','BRI','ABSI','AIP','HOMA_IR','CAR','MHR','TyG_hsCRP','hsCRP_BMI','hsCRP_WC','SII_BMI','NLR','PLR','SII','SIRI','PNI']
    rows.append({'record_id':r.index,'expression_id':r.index.removesuffix('-b'),
                 'base_id':r.component_terms,'arm_A_main':arm_a,'arm_B_24':arm_b,
                 'role':r.role,'exclusion_from_A':'' if arm_a else ('second GLM7 registration against full embedded base' if r.index=='GLM7' else 'log-scale sensitivity registration'),
                 'exclusion_from_B':'' if arm_b else 'outside the prespecified 24-expression demonstration',
                 'formula':r.formula,'expanded_formula':r.expanded_formula,
                 'contains_glucose':('glucose' in str(r.expanded_formula).lower() or 'TyG' in str(r.expanded_formula)),
                 'provenance':r.provenance,'endorsement':r.endorsement,
                 'audit_class_registered_base':r.reference_class})
manifest=pd.DataFrame(rows)
manifest.to_csv(out/'nhanes_record_manifest_v37.csv',index=False)
assert len(manifest)==26
assert manifest.expression_id.nunique()==25
assert manifest.arm_A_main.sum()==24
assert manifest.arm_B_24.sum()==24
assert manifest.arm_B_24.astype(bool).sum()==24
print(manifest[['record_id','arm_A_main','arm_B_24','contains_glucose']].to_string(index=False))
print('MANIFEST_OK: 26 records; 25 expressions; 24 Arm A; 24 Arm B')
