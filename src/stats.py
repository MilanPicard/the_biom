from scipy import stats as scistats
import numpy as np
def ttest(dfs,pvalues):
    stats = []
    groupbys = {}
    for i in range(len(dfs)):
        if i not in groupbys:
            groupby_i = dfs[i].groupby("box_category")
            groupbys[i] = groupby_i
        else:
            groupby_i = groupbys[i]

        groups_i = list(groupby_i.groups)
        for k in range(len(groups_i)):
            group_k = groupby_i.get_group(groups_i[k])
            for l in range(k+1,len(groups_i)):
                group_l = groupby_i.get_group(groups_i[l])
                r = scistats.ttest_ind(group_k["expression"],group_l["expression"])
                score = np.sum([1 if r.pvalue < p else 0 for p in pvalues])
                if(score>0):
                    stats.append((dfs[i].iloc[0]["gene"],groups_i[k],dfs[i].iloc[0]["gene"],groups_i[l],"*"*score))
        for j in range(i+1,len(dfs)):
            if j not in groupbys:
                groupby_j = dfs[j].groupby("box_category")
                groupbys[j] = groupby_j
            else:
                groupby_j = groupbys[j]
            groups_j = groupby_j.groups
            for k in range(len(groups_i)):
                group_k = groupby_i.get_group(groups_i[k])
                for group_l_id in groups_j.keys():#[groups_i[k]]:
                    group_l = groupby_j.get_group(group_l_id)
                    r = scistats.ttest_ind(group_k["expression"],group_l["expression"])
                    score = np.sum([1 if r.pvalue < p else 0 for p in pvalues])
                    if(score>0):
                        stats.append((dfs[i].iloc[0]["gene"],groups_i[k],dfs[j].iloc[0]["gene"],group_l_id,"*"*score))
    return stats