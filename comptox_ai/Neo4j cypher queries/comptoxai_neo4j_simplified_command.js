// Step 1: Identify all nodes linking PFNA, identified proteins, and DKD
MATCH (d:Disease {commonName: 'Adult Hepatocellular Carcinoma'})
MATCH (end:Gene)
WHERE end.geneSymbol IN ['PNPLA3']
CALL apoc.path.spanningTree(d, {
	relationshipFilter: '<GENEASSOCIATESWITHDISEASE|GENEINTERACTSWITHGENE',
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
WHERE end.geneSymbol IN ['PNPLA3']
CALL apoc.path.spanningTree(c, {
	relationshipFilter: 'CHEMICALINCREASESEXPRESSION>|CHEMICALDECREASESEXPRESSION>|GENEINTERACTSWITHGENE',
    minLevel: 1,
    maxLevel: 7, 
    endNodes: end
})
YIELD path
WITH nodes(path) as n, relationships(path) as p
RETURN n, p



// Step 2: Use those identified nodes and identify all relationships
MATCH (c:Chemical {xrefDTXSID: 'DTXSID8031863'})
MATCH (d:Disease {commonName: 'Diabetic Nephropathy'})
MATCH (node1:Gene) 
WHERE node1.geneSymbol IN ['ALPK1', 'APOA1', 'APOA2', 'APP', 'ATG3', 'CNDP1', 'CNDP2', 'COLGALT1', 'CREB3', 'FGF8', 'FGFR1', 'GDF11', 'GPNMB', 'HMGCS2', 'HSP90AA1', 'HSP90AB1', 'KPNA2', 'MSTN', 'NOS3', 'PON1', 'PSEN2', 'PSENEN', 'RNF114', 'SGK1', 'SMAD4', 'TGFB1', 'TGFBR1', 'TGFBR2', 'TXNIP', 'UBC'] 
WITH collect(id(node1))+collect(c)+collect(d) as nodes
CALL apoc.algo.cover(nodes)
YIELD rel
RETURN  startNode(rel), rel, endNode(rel);


// Step 2: Use those identified nodes and identify all relationships
MATCH (c:Chemical {xrefDTXSID: 'DTXSID8031863'})
MATCH (d:Disease {commonName: 'Diabetic Nephropathy'})
MATCH (node1:Gene) 
WHERE node1.geneSymbol IN ['ALPK1', 'APOA1', 'APOA2', 'APP', 'ATG3', 'CNDP1', 'CNDP2', 'COLGALT1', 'CREB3', 'FGF8', 'FGFR1', 'GDF11', 'GPNMB', 'HMGCS2', 'HSP90AA1', 'HSP90AB1', 'KPNA2', 'MSTN', 'NOS3', 'PON1', 'PSEN2', 'PSENEN', 'RNF114', 'SGK1', 'SMAD4', 'TGFB1', 'TGFBR1', 'TGFBR2', 'TXNIP', 'UBC'] 
WITH collect(id(node1))+collect(c)+collect(d) as nodes
CALL apoc.algo.cover(nodes)
YIELD rel
RETURN  startNode(rel), rel, endNode(rel);