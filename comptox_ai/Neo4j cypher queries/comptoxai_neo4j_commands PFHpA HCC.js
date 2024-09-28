

// HCC/PFAS single nodes alone:
MATCH (d:Disease {commonName: 'Adult Hepatocellular Carcinoma'}) RETURN d
MATCH (c:Chemical {xrefDTXSID: "DTXSID1037303"}) RETURN c
MATCH (g:Gene {geneSymbol: "COLGALT1"}) RETURN g

//Gene List
//'PNPLA3' // MEC
//'ACSL1', 'AKR1C1' // Ana's data

//eQTLs from open targets genetics
'PNPLA5', 'EFCAB6', 'RTL6', 'SULT4A1', 'PARVG', 
'SAMM50', 'MPPED1', 'PNPLA3', 'PARVB', 'SHISAL1'






// This returns a graph with all CHEMICALINCREASESEXPRESSION associations with PFNA
MATCH (c:Chemical)-[:CHEMICALINCREASESEXPRESSION|CHEMICALDECREASESEXPRESSION]-(Gene)
WHERE c.xrefDTXSID = 'DTXSID1037303' 
RETURN c, Gene

// Full command for identifying all nodes in the graph
MATCH (d:Disease {commonName: 'Adult Hepatocellular Carcinoma'})
MATCH (end:Gene)
WHERE end.geneSymbol IN ['PNPLA3', 'ACSL1', 'AKR1C1', 'PNPLA5', 'EFCAB6', 'RTL6', 'SULT4A1', 'PARVG', 'SAMM50', 'MPPED1', 'PARVB', 'SHISAL1']
CALL apoc.path.spanningTree(d, {
	relationshipFilter: '<GENEASSOCIATESWITHDISEASE|GENEINTERACTSWITHGENE|GENEINPATHWAY',
    minLevel: 1,
    maxLevel: 5, 
    endNodes: end
})
YIELD path
WITH nodes(path) as n, relationships(path) as p
RETURN n, p
UNION
MATCH (c:Chemical {xrefDTXSID: 'DTXSID1037303'})
MATCH (d:Disease {commonName: 'Adult Hepatocellular Carcinoma'})
MATCH (end:Gene)
WHERE end.geneSymbol IN ['PNPLA3', 'ACSL1', 'AKR1C1', 'PNPLA5', 'EFCAB6', 'RTL6', 'SULT4A1', 'PARVG', 'SAMM50', 'MPPED1', 'PARVB', 'SHISAL1']
CALL apoc.path.spanningTree(c, {
	relationshipFilter: 'CHEMICALINCREASESEXPRESSION>|CHEMICALDECREASESEXPRESSION>|GENEINTERACTSWITHGENE|GENEINPATHWAY',
    minLevel: 1,
    maxLevel: 7, 
    endNodes: end
})
YIELD path
WITH nodes(path) as n, relationships(path) as p
RETURN n, p






// Final query
MATCH (c:Chemical {xrefDTXSID: 'DTXSID8031863'}) 
MATCH (d:Disease {commonName: 'Liver carcinoma'})
MATCH (node1:Gene) WHERE node1.geneSymbol IN ['GSTA1', 'PGRMC1', 'RPS20', 'ACAA2', 'FGB', 'ANG', 'HERC5', 'FST', 'SKA3', 'TUBB2A', 'CENPU', 'UBD', 'PPARG', 'ANGPTL8', 'FAM72B', 'TP53', 'MAPK14', 'POR', 'FN1', 'SERPINE1', 'HADHA', 'PBK', 'ANGPTL4', 'HMGCS2', 'HADHB', 'SERPINA3', 'FABP1', 'SMARCD1', 'TK1', 'HAMP', 'UGT2B4', 'ARID1B', 'SUMO2', 'APOH', 'CYP3A7', 'JDP2', 'FBP1', 'CYP2B6', 'SLC25A47', 'ACAA1', 'PCLAF', 'EPS8L3', 'VKORC1', 'VCAM1', 'APOA2', 'UBC', 'STAT5A', 'AHSG', 'MT1X', 'CYP3A4', 'UGT2B7', 'FGA', 'ACTG1', 'LRRC59', 'APOA4', 'ACADVL', 'S100A6', 'PDK4', 'FATE1', 'TSLP', 'ACADM', 'AGT', 'SOCS3', 'MT1F', 'NFAT5', 'APOA5', 'DCN', 'UGT2B10', 'CYP2E1', 'CYP4A11', 'ALB', 'MMP9', 'COL4A1', 'F2', 'ADH1A', 'TPM1', 'AKR1B10', 'PLIN2', 'APP', 'TM4SF4', 'SAA4', 'GRB2', 'FGG', 'CYP2C9', 'ACSL1', 'PNPLA3', 'AKR1C3', 'ITGB1', 'CCL14', 'TIMP1', 'PTGR1', 'ACOX1', 'AKT2', 'ABCB4', 'RCHY1', 'APOF', 'CLEC4G', 'PCK1', 'RPL23', 'ESR1', 'CSNK2A2']
MATCH (node2:Pathway) WHERE node2.commonName IN ['Platelet activation, signaling and aggregation', 'Leukotriene metabolism', 'PPAR signaling pathway - Homo sapiens (human)', 'Regulation of lipid metabolism by Peroxisome proliferator-activated receptor alpha (PPARalpha)', 'Tyrosine metabolism', 'Aminosugars metabolism', 'Pathways in cancer - Homo sapiens (human)', 'Signal Transduction', 'Metabolism', 'Developmental Biology', 'Axon guidance', 'Metabolic pathways', 'Signaling by GPCR', 'Platelet degranulation ', 'Innate Immune System', 'PPARA activates gene expression', 'Fatty acid degradation - Homo sapiens (human)', 'Biological oxidations', 'Hemostasis', 'Extracellular matrix organization', 'Nuclear Receptors Meta-Pathway', 'Immune System', 'Proteoglycans in cancer - Homo sapiens (human)', 'Regulation of Insulin-like Growth Factor (IGF) transport and uptake by Insulin-like Growth Factor Binding Proteins (IGFBPs)', 'Metabolism of lipids and lipoproteins']
WITH collect(id(node1))+collect(c)+collect(d)+collect(id(node2)) as nodes 
CALL apoc.algo.cover(nodes) YIELD rel 
RETURN  startNode(rel), rel, endNode(rel);