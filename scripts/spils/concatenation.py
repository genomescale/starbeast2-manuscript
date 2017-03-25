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

def make_concat_xml(n_loci, data, run_name, log_every, chain_length, relaxed_clock = False, continuous_rates = False, species_tree_rates = False, adaptive_operators = False):
	log_every_str = str(log_every)
	chain_length_str = str(chain_length)

	gene_order = sorted(data)[:n_loci]
	first_gene = gene_order[0]
	species_set = set()

	for gene_symbol in gene_order:
		species_set.update(data[gene_symbol])

	n_extant_species = len(species_set)
	n_rates = (2 * n_extant_species) - 2
	n_rates_str = str(n_rates)
	one_on_n_str = str(1.0 / n_rates)

	species_order = sorted(species_set)

	beast = Element('beast', {'version': '2.0', 'namespace': 'beast.core:beast.evolution.alignment:beast.evolution.tree.coalescent:beast.core.util:beast.evolution.nuc:beast.evolution.operators:beast.evolution.sitemodel:beast.evolution.substitutionmodel:beast.evolution.likelihood'})
	beast_xml = xml.etree.ElementTree.ElementTree(beast)

	for gene_name in gene_order:
		gene_data = Element('data', {'id': gene_name, 'name': 'alignment'})
		for species_name in species_order:
			if species_name in data[gene_name]:
				sequence_name = gene_name + "_" + species_name
				sequence = data[gene_name][species_name]
				gene_species_tip_sequence = Element('sequence', {'totalcount': '4', 'id': sequence_name, 'value': sequence, 'taxon': species_name})
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
	taxonsuperset = Element('taxonset', {'id': 'taxonSuperSet', 'spec': 'TaxonSet', 'alignment': '@' + first_gene})
	species_tree.append(taxonsuperset)
	state.append(species_tree)

	for gene_name in gene_order:
		if relaxed_clock and not species_tree_rates:
			if continuous_rates:
				gene_branchrates_statenode = Element('parameter', {'id': 'branchRates:' + gene_name, 'name': 'stateNode', 'dimension': n_rates_str})
				gene_branchrates_statenode.text = one_on_n_str
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

	netdiversification_parameter = Element('parameter', {'upper': '1000.0', 'lower': '0.0', 'id': 'netDiversification', 'name': 'stateNode'})
	netdiversification_parameter.text = '100.0'
	state.append(netdiversification_parameter)
	mcmc_run.append(state)

	initializer = Element('init', {'spec': 'beast.util.TreeParser', 'initial': '@tree:species', 'taxa': '@' + first_gene, 'id': 'newick:species', 'IsLabelledNewick': 'true', 'newick': '((((s0:0.025,s1:0.025):0.025,s2:0.05):0.025,s3:0.075):0.025,s4:0.1);'})
	mcmc_run.append(initializer)

	posterior_distribution = Element('distribution', {'id': 'posterior', 'spec': 'util.CompoundDistribution'})
	prior_distribution = Element('distribution', {'id': 'prior', 'spec': 'util.CompoundDistribution'})
	species_likelihood_distribution = Element('distribution', {'birthDiffRate': '@netDiversification', 'spec': 'beast.evolution.speciation.BirthDeathGernhard08Model', 'tree': '@tree:species', 'id': 'birthDeathPrior', 'relativeDeathRate': '0.0'})
	netdiversificationprior = Element('prior', {'x': '@netDiversification', 'id': 'netDiversificationPrior', 'name': 'distribution'})
	netdiversificationprior_uniform = Element('Uniform', {'id': 'Uniform:netDiversificationPrior', 'lower': '0.0', 'upper': '1000.0', 'name': 'distr'})
	netdiversificationprior.append(netdiversificationprior_uniform)
	prior_distribution.append(species_likelihood_distribution)
	prior_distribution.append(netdiversificationprior)

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
		gene_likelihood_distribution = Element('distribution', {'data': '@' + gene_name, 'tree': '@tree:species', 'id': 'likelihood:' + gene_name, 'spec': 'TreeLikelihood'})
		gene_sitemodel = Element('siteModel', {'proportionInvariant': '0.0', 'mutationRate': '1.0', 'id': 'siteModel:' + gene_name, 'spec': 'SiteModel'})
		gene_jc_substmodel = Element('substModel', {'id': 'jc:' + gene_name, 'spec': 'JukesCantor'})
		gene_sitemodel.append(gene_jc_substmodel)
		gene_likelihood_distribution.append(gene_sitemodel)

		if not relaxed_clock:
			gene_clock_branchratemodel = Element('branchRateModel', {'clock.rate': '1.0', 'id': 'clock:' + gene_name, 'spec': 'beast.evolution.branchratemodel.StrictClockModel'})
		elif continuous_rates:
			branchratesprior = Element('prior', {'x': '@branchRates:' + gene_name, 'id': 'branchRatesPrior:' + gene_name, 'name': 'distribution'})
			branchratesprior_dirichlet = Element('Dirichlet', {'id': 'brancRatesPriorDistr:' + gene_name, 'name': 'distr'})
			branchrateconcentration = Element('parameter', {'id': 'branchRateConcentration:' + gene_name, 'name': 'alpha', 'estimate': 'false', 'dimension': n_rates_str})
			branchrateconcentration.text = '10.0'
			branchratesprior_dirichlet.append(branchrateconcentration)
			branchratesprior.append(branchratesprior_dirichlet)
			prior_distribution.append(branchratesprior)
			gene_clock_branchratemodel = Element('branchRateModel', {'clock.rate': '1.0', 'id': 'clock:' + gene_name, 'spec': 'starbeast2.DirichletRates', 'tree': '@tree:species', 'estimateRoot': 'false'})

			if species_tree_rates:
				gene_clock_branchratemodel.attrib['treeRates'] = '@branchRates:species'
			else:
				gene_clock_branchratemodel.attrib['treeRates'] = '@branchRates:' + gene_name
		elif species_tree_rates:
			gene_clock_branchratemodel = Element('branchRateModel', {'clock.rate': '1.0', 'id': 'clock:' + gene_name, 'spec': 'starbeast2.UncorrelatedRates', 'tree': '@tree:species', 'rates': '@branchRates:species', 'nBins': '100', 'estimateRoot': 'false', 'stdev': '0.16'})
		else:
			gene_clock_branchratemodel = Element('branchRateModel', {'clock.rate': '1.0', 'id': 'clock:' + gene_name, 'spec': 'starbeast2.UncorrelatedRates', 'tree': '@tree:species', 'nBins': '100', 'rates': '@branchRates:' + gene_name, 'estimateRoot': 'false', 'stdev': '0.32'})

		gene_likelihood_distribution.append(gene_clock_branchratemodel)
		likelihood_distribution.append(gene_likelihood_distribution)

	posterior_distribution.append(prior_distribution)
	posterior_distribution.append(likelihood_distribution)
	mcmc_run.append(posterior_distribution)

	netdiversificationscaler_operator = Element('operator', {'weight': '1.5', 'scaleFactor': '0.5', 'parameter': '@netDiversification', 'spec': 'ScaleOperator', 'id': 'netDiversificationScaler'})
	all_updown_operator = Element('operator', {'scaleFactor': '0.9', 'id': 'updown:all', 'weight': '25.0', 'spec': 'UpDownOperator'})
	netdiversification_up = Element('up', {'idref': 'netDiversification'})
	species_tree_down = Element('down', {'idref': 'tree:species'})
	all_updown_operator.append(netdiversification_up)
	all_updown_operator.append(species_tree_down)

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
					networkrateexchange_operator.attrib['weight']  = '5.0'
					diffuserateexchange_operator.attrib['weight']  = '5.0'
					childrateexchange_operator.attrib['weight']    = '5.0'
					childrenrateexchange_operator.attrib['weight'] = '5.0'
				else:
					networkrateexchange_operator.attrib['optimise']  = 'false'
					diffuserateexchange_operator.attrib['optimise']  = 'false'
					childrateexchange_operator.attrib['optimise']    = 'false'
					childrenrateexchange_operator.attrib['optimise'] = 'false'
					networkrateexchange_operator.attrib['weight']  = '15.0'
					diffuserateexchange_operator.attrib['weight']  = '15.0'
					childrateexchange_operator.attrib['weight']    = '15.0'
					childrenrateexchange_operator.attrib['weight'] = '15.0'
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
					branchrateswapper_operator.attrib['weight'] = '10.0'
					branchrateuniform_operator.attrib['weight'] = '10.0'
				else:
					branchrateswapper_operator.attrib['optimise'] = 'false'
					branchrateuniform_operator.attrib['optimise'] = 'false'
					branchrateswapper_operator.attrib['weight'] = '30.0'
					branchrateuniform_operator.attrib['weight'] = '30.0'
				mcmc_run.append(branchrateswapper_operator)
				mcmc_run.append(branchrateuniform_operator)

	species_treescaler_operator = Element('operator', {'tree': '@tree:species', 'scaleFactor': '0.99', 'id': 'treeScaler:species', 'weight': '10.0', 'spec': 'ScaleOperator'})
	species_treerootscaler_operator = Element('operator', {'rootOnly': 'true', 'scaleFactor': '0.9', 'tree': '@tree:species', 'spec': 'ScaleOperator', 'weight': '10.0', 'id': 'treeRootScaler:species'})
	species_uniform_operator = Element('operator', {'tree': '@tree:species', 'id': 'uniform:species', 'weight': '50.0', 'spec': 'Uniform'})
	mcmc_run.append(species_treescaler_operator)
	mcmc_run.append(species_treerootscaler_operator)
	mcmc_run.append(species_uniform_operator)

	for gene_name in gene_order:
		if relaxed_clock and not species_tree_rates:
			if continuous_rates:
				networkrateexchange_operator = Element('operator', {'id': 'networkRateExchange:' + gene_name, 'k': '2', 'treeRates': '@branchRates:' + gene_name, 'spec': 'starbeast2.NetworkRateExchange', 'delta': '0.5'})
				diffuserateexchange_operator = Element('operator', {'id': 'diffuseRateExchange:' + gene_name, 'k': '1', 'treeRates': '@branchRates:' + gene_name, 'spec': 'starbeast2.DiffuseRateExchange', 'delta': '0.5'})
				childrateexchange_operator = Element('operator', {'id': 'childRateExchange:' + gene_name, 'k': '1', 'treeRates': '@branchRates:' + gene_name, 'tree': '@tree:species', 'spec': 'starbeast2.ChildRateExchange', 'delta': '0.5'})
				childrenrateexchange_operator = Element('operator', {'id': 'childrenRateExchange:' + gene_name, 'k': '1', 'treeRates': '@branchRates:' + gene_name, 'tree': '@tree:species', 'spec': 'starbeast2.ChildrenRateExchange', 'delta': '0.5'})
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
	prior_log = Element('log', {'idref': 'prior'})
	birth_death_prior_log = Element('log', {'idref': 'birthDeathPrior'})
	likelihood_log = Element('log', {'idref': 'likelihood'})
	netdiversification_log = Element('log', {'idref': 'netDiversification'})
	speciestreeheight_log = Element('log', {'id': 'speciesTreeHeight', 'spec': 'beast.evolution.tree.TreeHeightLogger', 'tree': '@tree:species'})
	speciestreelength_log = Element('log', {'id': 'speciesTreeLength', 'spec': 'starbeast2.TreeLengthLogger', 'tree': '@tree:species'})
	branchrates_log = Element('log', {'idref': 'branchRates:species'})

	tracelog.append(posterior_log)
	tracelog.append(birth_death_prior_log)
	tracelog.append(prior_log)
	tracelog.append(likelihood_log)
	tracelog.append(netdiversification_log)
	tracelog.append(speciestreeheight_log)
	tracelog.append(speciestreelength_log)
	tracelog.append(branchrates_log)

	screenlog = Element('logger', {'model': '@posterior', 'logEvery': log_every_str, 'id': 'screenlog'})
	screenlog.append(posterior_log)
	screenlog.append(prior_log)
	screenlog.append(likelihood_log)

	species_treelog = Element('logger', {'fileName': run_name + '.trees', 'logEvery': log_every_str, 'mode': 'tree', 'id': 'treelog:species'})
	species_treelogger = Element('log', {'spec': 'beast.evolution.tree.TreeWithMetaDataLogger', 'id': 'treeLogger:species', 'tree': '@tree:species'})
	species_treelog.append(species_treelogger)

	if species_tree_rates:
		if continuous_rates:
			species_rate_model = Element('branchratemodel', {'clock.rate': '1.0', 'id': 'clock:species', 'spec': 'starbeast2.DirichletRates', 'tree': '@tree:species', 'estimateRoot': 'false', 'treeRates': '@branchRates:species', 'noCache': 'true'})
		else:
			species_rate_model = Element('branchratemodel', {'clock.rate': '1.0', 'id': 'clock:species', 'spec': 'starbeast2.UncorrelatedRates', 'tree': '@tree:species', 'estimateRoot': 'false', 'stdev': '0.16', 'nBins': '100', 'rates': '@branchRates:species', 'noCache': 'true'})
		species_treelogger.append(species_rate_model)

	mcmc_run.append(tracelog)
	mcmc_run.append(screenlog)
	mcmc_run.append(species_treelog)

	beast.append(mcmc_run)

	indent(beast)

	return beast_xml
