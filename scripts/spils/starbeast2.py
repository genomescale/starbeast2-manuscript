import collections
import xml.etree.ElementTree

Element = xml.etree.ElementTree.Element

def indent(elem, level = 0):
	i = '\n' + level * '  '
	if len(elem):
		if not elem.text or not elem.text.strip():
			elem.text = i + '  '
		if not elem.tail or not elem.tail.strip():
			elem.tail = i
		for elem in elem:
			indent(elem, level + 1)
		if not elem.tail or not elem.tail.strip():
			elem.tail = i
	else:
		if level and (not elem.tail or not elem.tail.strip()):
			elem.tail = i

def make_starbeast_xml(n_loci, data, run_name, log_every, chain_length, relaxed_clock = False, continuous_rates = False, species_tree_rates = False, analytical_pop_sizes = False, coordinated_topology = False, coordinated_heights = False, adaptive_operators = False):
	log_every_str = str(log_every)
	chain_length_str = str(chain_length)

	gene_order = sorted(data)[:n_loci]
	species_set = set()
	tip_sets = collections.defaultdict(set)

	tips_per_locus = collections.defaultdict(int)

	for gene_symbol in gene_order:
		species_data = data[gene_symbol]
		species_set.update(species_data)
		for species_name, tip_data in species_data.items():
			tips_per_locus[gene_symbol] += len(tip_data)
			tip_sets[species_name].update(tip_data)

	per_gene_rates = {}
	inverse_pgr = {}
	for gene_symbol, n_tips in tips_per_locus.items():
		n_gene_rates = 2 * (n_tips - 1)
		per_gene_rates[gene_symbol] = str(n_gene_rates)
		inverse_pgr[gene_symbol] = str(1.0 / n_gene_rates)

	n_extant_species = len(species_set)
	n_rates = (2 * n_extant_species) - 2
	n_rates_str = str(n_rates)
	one_on_n_str = str(1.0 / n_rates)

	species_order = sorted(species_set)
	tip_order = {}
	for species_name in species_order:
		tip_order[species_name] = sorted(tip_sets[species_name])

	beast = Element('beast', {'version': '2.0', 'namespace': 'beast.core:beast.evolution.alignment:beast.evolution.tree.coalescent:beast.core.util:beast.evolution.nuc:beast.evolution.operators:beast.evolution.sitemodel:beast.evolution.substitutionmodel:beast.evolution.likelihood'})
	beast_xml = xml.etree.ElementTree.ElementTree(beast)

	for gene_name in gene_order:
		gene_data = Element('data', {'id': gene_name, 'name': 'alignment'})
		for species_name in species_order:
			for tip_name in tip_order[species_name]:
				if tip_name in data[gene_name][species_name]:
					sequence_name = gene_name + "_" + tip_name
					sequence = data[gene_name][species_name][tip_name]
					gene_species_tip_sequence = Element('sequence', {'totalcount': '4', 'id': sequence_name, 'value': sequence, 'taxon': tip_name})
					gene_data.append(gene_species_tip_sequence)
		beast.append(gene_data)

	dirichlet_map = Element('map', {'name': 'Dirichlet'})
	uniform_map = Element('map', {'name': 'Uniform'})
	exponential_map = Element('map', {'name': 'Exponential'})
	lognormal_map = Element('map', {'name': 'LogNormal'})
	normal_map = Element('map', {'name': 'Normal'})
	beta_map = Element('map', {'name': 'Beta'})
	gamma_map = Element('map', {'name': 'Gamma'})
	laplacedistribution_map = Element('map', {'name': 'LaplaceDistribution'})
	prior_map = Element('map', {'name': 'prior'})
	inversegamma_map = Element('map', {'name': 'InverseGamma'})
	oneonx_map = Element('map', {'name': 'OneOnX'})
	dirichlet_map.text = 'beast.math.distributions.Dirichlet'
	uniform_map.text = 'beast.math.distributions.Uniform'
	exponential_map.text = 'beast.math.distributions.Exponential'
	lognormal_map.text = 'beast.math.distributions.LogNormalDistributionModel'
	normal_map.text = 'beast.math.distributions.Normal'
	beta_map.text = 'beast.math.distributions.Beta'
	gamma_map.text = 'beast.math.distributions.Gamma'
	laplacedistribution_map.text = 'beast.math.distributions.LaplaceDistribution'
	prior_map.text = 'beast.math.distributions.Prior'
	inversegamma_map.text = 'beast.math.distributions.InverseGamma'
	oneonx_map.text = 'beast.math.distributions.OneOnX'
	beast.append(dirichlet_map)
	beast.append(uniform_map)
	beast.append(exponential_map)
	beast.append(lognormal_map)
	beast.append(normal_map)
	beast.append(beta_map)
	beast.append(gamma_map)
	beast.append(laplacedistribution_map)
	beast.append(prior_map)
	beast.append(inversegamma_map)
	beast.append(oneonx_map)

	mcmc_run = Element('run', {'chainLength': chain_length_str, 'storeEvery': log_every_str, 'id': 'mcmc', 'spec': 'MCMC'})
	state = Element('state', {'storeEvery': log_every_str, 'id': 'state'})
	species_tree = Element('tree', {'id': 'tree:species', 'name': 'stateNode'})
	taxonsuperset = Element('taxonset', {'id': 'taxonSuperSet', 'spec': 'TaxonSet'})
	for species_name in species_order:
		species_taxon = Element('taxon', {'id': species_name, 'spec': 'TaxonSet'})
		for tip_name in tip_order[species_name]:
			species_tip_taxon = Element('taxon', {'id': tip_name, 'spec': 'Taxon'})
			species_taxon.append(species_tip_taxon)
		taxonsuperset.append(species_taxon)
	species_tree.append(taxonsuperset)
	state.append(species_tree)

	for gene_name in gene_order:
		gene_tree = Element('tree', {'id': 'tree:' + gene_name, 'name': 'stateNode'})
		gene_taxonset = Element('taxonset', {'id': 'TaxonSet:' + gene_name, 'alignment': '@' + gene_name, 'spec': 'TaxonSet'})
		gene_tree.append(gene_taxonset)
		state.append(gene_tree)

		if relaxed_clock and not species_tree_rates:
			if continuous_rates:
				gene_branchrates_statenode = Element('parameter', {'id': 'branchRates:' + gene_name, 'name': 'stateNode', 'dimension': per_gene_rates[gene_name]})
				gene_branchrates_statenode.text = inverse_pgr[gene_name]
				state.append(gene_branchrates_statenode)
			else:
				gene_branchrates_statenode = Element('stateNode', {'id': 'branchRates:' + gene_name, 'spec': 'parameter.IntegerParameter', 'dimension': '100'})
				gene_branchrates_statenode.text = '50'
				state.append(gene_branchrates_statenode)

	if relaxed_clock and species_tree_rates:
		if continuous_rates:
			species_branchrates_parameter = Element('parameter', {'id': 'branchRates:species', 'name': 'stateNode', 'lower': '0.0', 'dimension': n_rates_str})
			species_branchrates_parameter.text = one_on_n_str
			state.append(species_branchrates_parameter)
		else:
			species_branchrates_statenode = Element('stateNode', {'id': 'branchRates:species', 'spec': 'parameter.IntegerParameter', 'dimension': '100'})
			species_branchrates_statenode.text = '50'
			state.append(species_branchrates_statenode)

	if not analytical_pop_sizes:
		popsizes_parameter = Element('parameter', {'lower': '0.0', 'id': 'popSizes', 'name': 'stateNode'})
		popsizes_parameter.text = '0.001'
		state.append(popsizes_parameter)

	netdiversification_parameter = Element('parameter', {'upper': '1000.0', 'lower': '0.0', 'id': 'netDiversification', 'name': 'stateNode'})
	popmean_parameter = Element('parameter', {'lower': '0.0', 'id': 'popMean', 'name': 'stateNode'})
	netdiversification_parameter.text = '100.0'
	popmean_parameter.text = '1.0'
	state.append(netdiversification_parameter)
	state.append(popmean_parameter)
	mcmc_run.append(state)

	initializer = Element('init', {'spec': 'starbeast2.StarBeastInitializer', 'estimate': 'false', 'method': 'random', 'speciesTree': '@tree:species', 'id': 'initializer', 'newick': '((((s0:0.025,s1:0.025):0.025,s2:0.05):0.025,s3:0.075):0.025,s4:0.1);'})
	for gene_name in gene_order:
		gene_tree = Element('geneTree', {'idref': 'tree:' + gene_name})
		initializer.append(gene_tree)
	if analytical_pop_sizes:
		popmodel = Element('populationModel', {'spec': 'starbeast2.ConstantPopulationIO', 'populationShape': '3.0', 'populationMean': '@popMean', 'id': 'popModel'})
	else:
		popmodel = Element('populationModel', {'spec': 'starbeast2.ConstantPopulation', 'populationSizes': '@popSizes', 'id': 'popModel'})
	initializer.append(popmodel)
	mcmc_run.append(initializer)

	posterior_distribution = Element('distribution', {'id': 'posterior', 'spec': 'util.CompoundDistribution'})
	speciescoalescent_distribution = Element('distribution', {'populationModel': '@popModel', 'id': 'speciescoalescent', 'spec': 'starbeast2.MultispeciesCoalescent'})
	speciestree = Element('speciesTree', {'tree': '@tree:species', 'id': 'speciesTree', 'spec': 'starbeast2.SpeciesTree'})
	speciescoalescent_distribution.append(speciestree)
	for gene_name in gene_order:
		gene_genetree = Element('geneTree', {'spec': 'starbeast2.GeneTree', 'ploidy': '2.0', 'tree': '@tree:' + gene_name, 'id': 'geneTree:' + gene_name, 'speciesTree': '@speciesTree'})
		speciescoalescent_distribution.append(gene_genetree)

	prior_distribution = Element('distribution', {'id': 'prior', 'spec': 'util.CompoundDistribution'})
	species_likelihood_distribution = Element('distribution', {'birthDiffRate': '@netDiversification', 'spec': 'beast.evolution.speciation.BirthDeathGernhard08Model', 'tree': '@tree:species', 'id': 'birthDeathPrior', 'relativeDeathRate': '0.0'})
	popmeanprior = Element('prior', {'x': '@popMean', 'id': 'popMeanPrior', 'name': 'distribution'})
	netdiversificationprior = Element('prior', {'x': '@netDiversification', 'id': 'netDiversificationPrior', 'name': 'distribution'})
	popmeanprior_oneonx = Element('OneOnX', {'id': 'OneOnX:popMeanPrior', 'name': 'distr'})
	netdiversificationprior_uniform = Element('Uniform', {'id': 'Uniform:netDiversificationPrior', 'lower': '0.0', 'upper': '1000.0', 'name': 'distr'})
	popmeanprior.append(popmeanprior_oneonx)
	netdiversificationprior.append(netdiversificationprior_uniform)
	prior_distribution.append(species_likelihood_distribution)
	prior_distribution.append(popmeanprior)
	prior_distribution.append(netdiversificationprior)

	if not analytical_pop_sizes:
		popsizesprior = Element('prior', {'x': '@popSizes', 'id': 'popSizesPrior', 'name': 'distribution'})
		popsizesprior_inversegamma = Element('distr', {'alpha': '3.0', 'mean': '@popMean', 'id': 'AltInverseGamma:popSizesPrior', 'spec': 'starbeast2.AltInverseGamma'})
		popsizesprior.append(popsizesprior_inversegamma)
		prior_distribution.append(popsizesprior)

	if relaxed_clock and species_tree_rates and continuous_rates:
		branchratesprior = Element('prior', {'x': '@branchRates:species', 'id': 'branchRatesPrior', 'name': 'distribution'})
		branchratesprior_dirichlet = Element('Dirichlet', {'id': 'Dirichlet:branchRates', 'name': 'distr'})
		branchrateconcentration = Element('parameter', {'id': 'branchRateConcentration', 'name': 'alpha', 'estimate': 'false', 'dimension': n_rates_str})
		branchrateconcentration.text = '40.0'
		branchratesprior_dirichlet.append(branchrateconcentration)
		branchratesprior.append(branchratesprior_dirichlet)
		prior_distribution.append(branchratesprior)

	likelihood_distribution = Element('distribution', {'id': 'likelihood', 'spec': 'util.CompoundDistribution'})
	for gene_i, gene_name in enumerate(gene_order):
		gene_likelihood_distribution = Element('distribution', {'data': '@' + gene_name, 'tree': '@tree:' + gene_name, 'id': 'likelihood:' + gene_name, 'spec': 'TreeLikelihood'})
		gene_sitemodel = Element('siteModel', {'proportionInvariant': '0.0', 'mutationRate': '1.0', 'id': 'siteModel:' + gene_name, 'spec': 'SiteModel'})
		gene_jc_substmodel = Element('substModel', {'id': 'jc:' + gene_name, 'spec': 'JukesCantor'})
		gene_sitemodel.append(gene_jc_substmodel)
		gene_likelihood_distribution.append(gene_sitemodel)

		if not relaxed_clock:
			gene_clock_branchratemodel = Element('branchRateModel', {'clock.rate': '1.0', 'id': 'clock:' + gene_name, 'spec': 'beast.evolution.branchratemodel.StrictClockModel'})

		elif species_tree_rates:
			if gene_i > 0:
				gene_clock_branchratemodel = Element('branchRateModel', {'clock.rate': '1.0', 'id': 'clock:' + gene_name, 'spec': 'starbeast2.StarBeastClock', 'geneTree': '@geneTree:' + gene_name, 'speciesTreeRates': '@clock:species'})
			elif continuous_rates:
				gene_clock_branchratemodel = Element('branchRateModel', {'clock.rate': '1.0', 'id': 'clock:' + gene_name, 'spec': 'starbeast2.StarBeastClock', 'geneTree': '@geneTree:' + gene_name})
				speciestreerates = Element('speciesTreeRates', {'id': 'clock:species', 'treeRates': '@branchRates:species', 'spec': 'starbeast2.DirichletRates', 'tree': '@tree:species', 'estimateRoot': 'false'})
				gene_clock_branchratemodel.append(speciestreerates)
			else:
				gene_clock_branchratemodel = Element('branchRateModel', {'clock.rate': '1.0', 'id': 'clock:' + gene_name, 'spec': 'starbeast2.StarBeastClock', 'geneTree': '@geneTree:' + gene_name})
				speciestreerates = Element('speciesTreeRates', {'id': 'clock:species', 'nBins': '100', 'rates': '@branchRates:species', 'spec': 'starbeast2.UncorrelatedRates', 'tree': '@tree:species', 'estimateRoot': 'false', 'stdev': '0.16'})
				gene_clock_branchratemodel.append(speciestreerates)

		elif continuous_rates:
			branchratesprior = Element('prior', {'x': '@branchRates:' + gene_name, 'id': 'branchRatesPrior:' + gene_name, 'name': 'distribution'})
			branchratesprior_dirichlet = Element('Dirichlet', {'id': 'brancRatesPriorDistr:' + gene_name, 'name': 'distr'})
			branchrateconcentration = Element('parameter', {'id': 'branchRateConcentration:' + gene_name, 'name': 'alpha', 'estimate': 'false', 'dimension': per_gene_rates[gene_name]})
			branchrateconcentration.text = '10.0'
			branchratesprior_dirichlet.append(branchrateconcentration)
			branchratesprior.append(branchratesprior_dirichlet)
			prior_distribution.append(branchratesprior)
			gene_clock_branchratemodel = Element('branchRateModel', {'clock.rate': '1.0', 'id': 'clock:' + gene_name, 'spec': 'starbeast2.DirichletRates', 'tree': '@tree:' + gene_name, 'treeRates': '@branchRates:' + gene_name, 'estimateRoot': 'false'})
		else:
			gene_clock_branchratemodel = Element('branchRateModel', {'clock.rate': '1.0', 'id': 'clock:' + gene_name, 'spec': 'starbeast2.UncorrelatedRates', 'tree': '@tree:' + gene_name, 'nBins': '100', 'rates': '@branchRates:' + gene_name, 'estimateRoot': 'false', 'stdev': '0.32'})

		gene_likelihood_distribution.append(gene_clock_branchratemodel)
		likelihood_distribution.append(gene_likelihood_distribution)

	posterior_distribution.append(speciescoalescent_distribution)
	posterior_distribution.append(prior_distribution)
	posterior_distribution.append(likelihood_distribution)
	mcmc_run.append(posterior_distribution)

	popmeanscaler_operator = Element('operator', {'parameter': '@popMean', 'scaleFactor': '0.75', 'id': 'popMeanScaler', 'weight': '12.0', 'spec': 'ScaleOperator'})
	netdiversificationscaler_operator = Element('operator', {'weight': '6.0', 'scaleFactor': '0.5', 'parameter': '@netDiversification', 'spec': 'ScaleOperator', 'id': 'netDiversificationScaler'})
	all_updown_operator = Element('operator', {'scaleFactor': '0.9', 'id': 'updown:all', 'weight': '100.0', 'spec': 'UpDownOperator'})
	netdiversification_up = Element('up', {'idref': 'netDiversification'})
	popmean_down = Element('down', {'idref': 'popMean'})
	species_tree_down = Element('down', {'idref': 'tree:species'})
	all_updown_operator.append(netdiversification_up)
	all_updown_operator.append(popmean_down)
	all_updown_operator.append(species_tree_down)
	for gene_name in gene_order:
		gene_tree_down = Element('down', {'idref': 'tree:' + gene_name})
		all_updown_operator.append(gene_tree_down)

	if not analytical_pop_sizes:
		popsizescaler_operator = Element('operator', {'parameter': '@popSizes', 'scaleFactor': '0.5', 'id': 'popSizeScaler', 'weight': '40.0', 'spec': 'ScaleOperator'})
		popsizeswapper_operator = Element('operator', {'parameter': '@popSizes', 'id': 'popSizeSwapper', 'weight': '40.0', 'spec': 'SwapOperator'})
		popsizes_down = Element('down', {'idref': 'popSizes'})
		all_updown_operator.append(popsizes_down)
		mcmc_run.append(popsizescaler_operator)
		mcmc_run.append(popsizeswapper_operator)

	mcmc_run.append(popmeanscaler_operator)
	mcmc_run.append(netdiversificationscaler_operator)
	mcmc_run.append(all_updown_operator)

	if relaxed_clock:
		if species_tree_rates:
			if continuous_rates:
				networkrateexchange_operator = Element('operator', {'id': 'networkRateExchange:species', 'k': '2', 'treeRates': '@branchRates:species', 'spec': 'starbeast2.NetworkRateExchange', 'delta': '0.5'})
				diffuserateexchange_operator = Element('operator', {'id': 'diffuseRateExchange:species', 'k': '1', 'treeRates': '@branchRates:species', 'spec': 'starbeast2.DiffuseRateExchange', 'delta': '0.5'})
				childrateexchange_operator = Element('operator', {'id': 'childRateExchange:species', 'k': '1', 'treeRates': '@branchRates:species', 'tree': '@tree:species', 'spec': 'starbeast2.ChildRateExchange', 'delta': '0.5'})
				childrenrateexchange_operator = Element('operator', {'id': 'childrenRateExchange:species', 'k': '1', 'treeRates': '@branchRates:species', 'tree': '@tree:species', 'spec': 'starbeast2.ChildrenRateExchange', 'delta': '0.5'})
				if adaptive_operators:
					networkrateexchange_operator.attrib['optimise']  = 'true'
					diffuserateexchange_operator.attrib['optimise']  = 'true'
					childrateexchange_operator.attrib['optimise']    = 'true'
					childrenrateexchange_operator.attrib['optimise'] = 'true'
					networkrateexchange_operator.attrib['weight']  = '20.0'
					diffuserateexchange_operator.attrib['weight']  = '20.0'
					childrateexchange_operator.attrib['weight']    = '20.0'
					childrenrateexchange_operator.attrib['weight'] = '20.0'
				else:
					networkrateexchange_operator.attrib['optimise']  = 'false'
					diffuserateexchange_operator.attrib['optimise']  = 'false'
					childrateexchange_operator.attrib['optimise']    = 'false'
					childrenrateexchange_operator.attrib['optimise'] = 'false'
					networkrateexchange_operator.attrib['weight']  = '60.0'
					diffuserateexchange_operator.attrib['weight']  = '60.0'
					childrateexchange_operator.attrib['weight']    = '60.0'
					childrenrateexchange_operator.attrib['weight'] = '60.0'
				mcmc_run.append(networkrateexchange_operator)
				mcmc_run.append(diffuserateexchange_operator)
				mcmc_run.append(childrateexchange_operator)
				mcmc_run.append(childrenrateexchange_operator)
			else:
				branchrateswapper_operator = Element('operator', {'id': 'branchRateCycle:species', 'k': '2', 'treeRates': '@branchRates:species', 'spec': 'starbeast2.DiscreteRateCycle'})
				branchrateuniform_operator = Element('operator', {'id': 'branchRateUniform:species', 'k': '1', 'treeRates': '@branchRates:species', 'spec': 'starbeast2.DiscreteRateUniform'})
				if adaptive_operators:
					branchrateswapper_operator.attrib['optimise'] = 'true'
					branchrateuniform_operator.attrib['optimise'] = 'true'
					branchrateswapper_operator.attrib['weight'] = '40.0'
					branchrateuniform_operator.attrib['weight'] = '40.0'
				else:
					branchrateswapper_operator.attrib['optimise'] = 'false'
					branchrateuniform_operator.attrib['optimise'] = 'false'
					branchrateswapper_operator.attrib['weight'] = '120.0'
					branchrateuniform_operator.attrib['weight'] = '120.0'
				mcmc_run.append(branchrateswapper_operator)
				mcmc_run.append(branchrateuniform_operator)

	species_treescaler_operator = Element('operator', {'tree': '@tree:species', 'scaleFactor': '0.99', 'id': 'treeScaler:species', 'weight': '40.0', 'spec': 'ScaleOperator'})
	species_treerootscaler_operator = Element('operator', {'rootOnly': 'true', 'scaleFactor': '0.9', 'tree': '@tree:species', 'spec': 'ScaleOperator', 'weight': '40.0', 'id': 'treeRootScaler:species'})
	species_uniform_operator = Element('operator', {'tree': '@tree:species', 'id': 'uniform:species', 'weight': '200.0', 'spec': 'Uniform'})
	mcmc_run.append(species_treescaler_operator)
	mcmc_run.append(species_treerootscaler_operator)
	mcmc_run.append(species_uniform_operator)

	if coordinated_heights:
		coordinated_heights_operator = Element('operator', {'tree': '@tree:species', 'speciesTree': '@speciesTree', 'id': 'coordinatedUniform', 'weight': '200.0', 'spec': 'starbeast2.CoordinatedUniform'})
		coordinated_exponential_operator = Element('operator', {'tree': '@tree:species', 'speciesTree': '@speciesTree', 'id': 'coordinatedExponential', 'weight': '200.0', 'spec': 'starbeast2.CoordinatedExponential', 'beta': '0.003'})
		mcmc_run.append(coordinated_heights_operator)
		mcmc_run.append(coordinated_exponential_operator)

	for gene_name in gene_order:
		coordinated_operator_genetree = Element('geneTree', {'idref': 'tree:' + gene_name})
		if coordinated_topology:
			species_wide_operator.append(coordinated_operator_genetree)
			species_narrow_operator.append(coordinated_operator_genetree)
		if coordinated_heights:
			coordinated_heights_operator.append(coordinated_operator_genetree)
			coordinated_exponential_operator.append(coordinated_operator_genetree)

		gene_treescaler_operator = Element('operator', {'tree': '@tree:' + gene_name, 'scaleFactor': '0.95', 'id': 'treeScaler:' + gene_name, 'weight': '3.0', 'spec': 'ScaleOperator'})
		gene_treerootscaler_operator = Element('operator', {'weight': '3.0', 'rootOnly': 'true', 'scaleFactor': '0.7', 'tree': '@tree:' + gene_name, 'spec': 'ScaleOperator', 'id': 'treeRootScaler:' + gene_name})
		gene_uniform_operator = Element('operator', {'tree': '@tree:' + gene_name, 'id': 'uniform:' + gene_name, 'weight': '15.0', 'spec': 'Uniform'})
		gene_subtreeslide_operator = Element('operator', {'spec': 'SubtreeSlide', 'tree': '@tree:' + gene_name, 'id': 'subtreeSlide:' + gene_name, 'weight': '15.0', 'size': '0.002'})
		gene_narrow_operator = Element('operator', {'tree': '@tree:' + gene_name, 'id': 'narrow:' + gene_name, 'weight': '15.0', 'spec': 'Exchange'})
		gene_wide_operator = Element('operator', {'spec': 'Exchange', 'tree': '@tree:' + gene_name, 'isNarrow': 'false', 'weight': '15.0', 'id': 'wide:' + gene_name})
		gene_wilsonbalding_operator = Element('operator', {'tree': '@tree:' + gene_name, 'id': 'WilsonBalding:' + gene_name, 'weight': '15.0', 'spec': 'WilsonBalding'})
		mcmc_run.append(gene_treescaler_operator)
		mcmc_run.append(gene_treerootscaler_operator)
		mcmc_run.append(gene_uniform_operator)
		mcmc_run.append(gene_subtreeslide_operator)
		mcmc_run.append(gene_narrow_operator)
		mcmc_run.append(gene_wide_operator)
		mcmc_run.append(gene_wilsonbalding_operator)

		if relaxed_clock and not species_tree_rates:
			if continuous_rates:
				networkrateexchange_operator = Element('operator', {'id': 'networkRateExchange:' + gene_name, 'k': '2', 'treeRates': '@branchRates:' + gene_name, 'spec': 'starbeast2.NetworkRateExchange', 'delta': '0.5'})
				diffuserateexchange_operator = Element('operator', {'id': 'diffuseRateExchange:' + gene_name, 'k': '1', 'treeRates': '@branchRates:' + gene_name, 'spec': 'starbeast2.DiffuseRateExchange', 'delta': '0.5'})
				childrateexchange_operator = Element('operator', {'id': 'childRateExchange:' + gene_name, 'k': '1', 'treeRates': '@branchRates:' + gene_name, 'tree': '@tree:' + gene_name, 'spec': 'starbeast2.ChildRateExchange', 'delta': '0.5'})
				childrenrateexchange_operator = Element('operator', {'id': 'childrenRateExchange:' + gene_name, 'k': '1', 'treeRates': '@branchRates:' + gene_name, 'tree': '@tree:' + gene_name, 'spec': 'starbeast2.ChildrenRateExchange', 'delta': '0.5'})
				if adaptive_operators:
					networkrateexchange_operator.attrib['optimise']  = 'true'
					diffuserateexchange_operator.attrib['optimise']  = 'true'
					childrateexchange_operator.attrib['optimise']    = 'true'
					childrenrateexchange_operator.attrib['optimise'] = 'true'
					networkrateexchange_operator.attrib['weight']  = '1.5'
					diffuserateexchange_operator.attrib['weight']  = '1.5'
					childrateexchange_operator.attrib['weight']    = '1.5'
					childrenrateexchange_operator.attrib['weight'] = '1.5'
				else:
					networkrateexchange_operator.attrib['optimise']  = 'false'
					diffuserateexchange_operator.attrib['optimise']  = 'false'
					childrateexchange_operator.attrib['optimise']    = 'false'
					childrenrateexchange_operator.attrib['optimise'] = 'false'
					networkrateexchange_operator.attrib['weight']  = '4.5'
					diffuserateexchange_operator.attrib['weight']  = '4.5'
					childrateexchange_operator.attrib['weight']    = '4.5'
					childrenrateexchange_operator.attrib['weight'] = '4.5'
				mcmc_run.append(networkrateexchange_operator)
				mcmc_run.append(diffuserateexchange_operator)
				mcmc_run.append(childrateexchange_operator)
				mcmc_run.append(childrenrateexchange_operator)
			else:
				branchrateswapper_operator = Element('operator', {'id': 'branchRateCycle:' + gene_name, 'k': '2', 'treeRates': '@branchRates:' + gene_name, 'spec': 'starbeast2.DiscreteRateCycle'})
				branchrateuniform_operator = Element('operator', {'id': 'branchRateUniform:' + gene_name, 'k': '1', 'treeRates': '@branchRates:' + gene_name, 'spec': 'starbeast2.DiscreteRateUniform'})
				if adaptive_operators:
					branchrateswapper_operator.attrib['optimise'] = 'true'
					branchrateuniform_operator.attrib['optimise'] = 'true'
					branchrateswapper_operator.attrib['weight'] = '3.0'
					branchrateuniform_operator.attrib['weight'] = '3.0'
				else:
					branchrateswapper_operator.attrib['optimise'] = 'false'
					branchrateuniform_operator.attrib['optimise'] = 'false'
					branchrateswapper_operator.attrib['weight'] = '9.0'
					branchrateuniform_operator.attrib['weight'] = '9.0'
				mcmc_run.append(branchrateswapper_operator)
				mcmc_run.append(branchrateuniform_operator)

	tracelog = Element('logger', {'sort': 'smart', 'model': '@posterior', 'logEvery': log_every_str, 'id': 'tracelog', 'fileName': run_name + '.log'})
	posterior_log = Element('log', {'idref': 'posterior'})
	speciescoalescent_log = Element('log', {'idref': 'speciescoalescent'})
	prior_log = Element('log', {'idref': 'prior'})
	birth_death_prior_log = Element('log', {'idref': 'birthDeathPrior'})
	likelihood_log = Element('log', {'idref': 'likelihood'})
	netdiversification_log = Element('log', {'idref': 'netDiversification'})
	popmean_log = Element('log', {'idref': 'popMean'})
	speciestreeheight_log = Element('log', {'id': 'speciesTreeHeight', 'spec': 'beast.evolution.tree.TreeHeightLogger', 'tree': '@tree:species'})
	speciestreelength_log = Element('log', {'id': 'speciesTreeLength', 'spec': 'starbeast2.TreeLengthLogger', 'tree': '@tree:species'})
	branchrates_log = Element('log', {'idref': 'branchRates:species'})

	tracelog.append(posterior_log)
	tracelog.append(speciescoalescent_log)
	tracelog.append(birth_death_prior_log)
	tracelog.append(prior_log)
	tracelog.append(likelihood_log)
	tracelog.append(netdiversification_log)
	tracelog.append(popmean_log)
	tracelog.append(speciestreeheight_log)
	tracelog.append(speciestreelength_log)
	tracelog.append(branchrates_log)

	screenlog = Element('logger', {'model': '@posterior', 'logEvery': log_every_str, 'id': 'screenlog'})
	screenlog.append(posterior_log)
	screenlog.append(speciescoalescent_log)
	screenlog.append(prior_log)
	screenlog.append(likelihood_log)

	species_treelog = Element('logger', {'fileName': run_name + '.trees', 'logEvery': log_every_str, 'mode': 'tree', 'id': 'treelog:species'})
	species_treelogger = Element('log', {'spec': 'starbeast2.SpeciesTreeLogger', 'id': 'treeLogger:species', 'speciesTree': '@speciesTree'})
	species_treelog.append(species_treelogger)

	if not analytical_pop_sizes:
		species_treelogger.attrib['populationmodel'] = '@popModel'

	if species_tree_rates:
		species_treelogger.attrib['branchratemodel'] = '@clock:species'

	mcmc_run.append(tracelog)
	mcmc_run.append(screenlog)
	mcmc_run.append(species_treelog)

	# for gene_name in gene_order:
	# 	output_filename = run_name + '.' + gene_name + '.trees'
	# 	gene_treelog = Element('logger', {'fileName': output_filename, 'logEvery': log_every_str, 'mode': 'tree', 'id': 'treelog:' + gene_name})
	# 	gene_treelogger = Element('log', {'tree': '@tree:' + gene_name, 'id': 'treeLogger:' + gene_name, 'spec': 'beast.evolution.tree.TreeWithMetaDataLogger'})
	# 	gene_treelog.append(gene_treelogger)

	# 	if relaxed_clock:
	# 		gene_treelogger.attrib['branchratemodel'] = '@clock:' + gene_name

	# 	mcmc_run.append(gene_treelog)

	beast.append(mcmc_run)

	indent(beast)

	return beast_xml
