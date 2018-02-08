# mutations

Deepvariant相关：  
mutsites_dv：从bam文件中提取候选突变位点  
addRef_dv：给候选突变位点加上图片覆盖的ref，顺便提供bam文件中部分位点缺失的del区域ref序列  
encode_dv：生成候选突变位点的图像数据  
train_dv：训练以及验证  
  
其他：  
aInb：判断两数组关系  
algorithm：根据突变位点特征判断突变  
refReader：按行获取给定候选位点的对应图像覆盖的ref  
refReader2：按染色体获取给定候选位点的对应图像覆盖的ref  
