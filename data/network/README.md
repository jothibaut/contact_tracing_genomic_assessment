For a matter of data privacy, those file do not reflect the real contact tracing network used in this article.
Those files show an example of how we constructed the contact tracing network.
The 'pairs' files represent social links between the patients:
	- source_id, target_id: ID's of patients within the contact tracing network.
	- EPI.x, EPI.y, variant.x,  variant.y, strain.x, strain.y: Genomic data linked to source (.x) and target patients (.y), as available within GISAID.
The 'component' files represent the contact groups within this network.
	- node: the strain id, as stored within 'strain.x' and 'strain.y'
	- component: the contact group, within the contact tracing network. Strains that are interconnected within the network should have the same component number.