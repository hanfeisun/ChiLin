from chilin.dc import *

jinja_env = Environment(loader = PackageLoader('chilin', 'template'),
                        block_start_string = '\BLOCK{',
                        block_end_string = '}',
                        variable_start_string = '\VAR{',
                        variable_end_string = '}',
                        comment_start_string = '\#{',
                        comment_end_string = '}',
                        line_statement_prefix = '%-',
                        line_comment_prefix = '%#',
                        trim_blocks = True,
                        autoescape = False,
                        )
class QC_Controller(object):
	"""
	All the class in the module derives from this class
	"""
	def __init__(self):
#		self.pathfinder = PathFinder(conf)
		self.env = jinja_env
		print 'pass'
		self.template = self.env.get_template('template.tex')
		self.has_run = False
		print 'QC control'

	def run(self):
		""" Run some QC tools or do some time-costing statistics """
		self.has_run = True
		return True

	def _check(self):
		""" Check whether the quality of the dataset is ok. """
		if not self.has_run:
			self.run()
		return True

	def _render(self, template = None):
		""" Generate the latex code for current section. """
#		self.template.render({})
		pass


class RawQC(QC_Controller):
	"""  
	RawQC aims to perform some simple quality control checks to ensure that the raw data looks good and there are no problems or biases in your data.
	"""
	def __init__(self):
		super(RawQC, self).__init__()
		print 'init basic_qc'
	def _basic_info(self):
		self.basic_stat = ['basic_info']
		print self.basic_stat
		""" Basic description of the ChIP-seq  input dataset """
		return self.basic_stat

	def _fastqc_info(self):
		self.has_fastqc = True
		""" QC analysis of the raw Chip-seq data, including sequence quality score of particularity raw data and the cumulative percentage plot of the sequence quality scores of all historic data.
		"""
		print 'fastqc\n'          
		return ['fastqc_summary'],'fastqc_pdf'
	def run(self):
		""" Run some RawQC functions to get final result."""
		self.RawQC_check = True
		self.basic_stat = self._basic_info() #list format
		self.fastqc_summary_stat,self.fastqc_graph_stat = self._fastqc_info()
		self._check()
		self._render()
	def _check(self):
		"""
		Check whether the FastQC's result is ok
		"""
		if self.has_fastqc:
			print 'input self.fastqc_summary_stat for judge'
		else:
			print ' no need fastqc'
	def _render(self):
		print self.fastqc_summary_stat
		temp = self.template.render(RawQC_check = self.RawQC_check,prefix_datasetid = 'id',basic_table = self.basic_stat,fastqc_check = self.has_fastqc,fastqc_graph = self.fastqc_graph_stat)
		print temp 


		
class MappingQC(QC_Controller):
	""" MappingQC aims to describe the mapping quality of the sequence alignment. """
	def __init__(self):
		super(MappingQC, self).__init__()		
		print 'init mapping qc'

	def _basic_mapping_statistics_info(self):
		""" Stastic summary of mapping result for each sample. """
		print 'basic_mapping_statistics'
		self.mappable_summary_stat = 'basic_mapping_table'
		return self.mappable_summary_stat

		""" Cumulative percentage plot to  describe the  mappable ratio quality of all historic data. """
	def _mappable_ratio_info(self):
		""" Cumulative percentage plot to  describe the  mappable ratio quality of all historic data."""
		print 'mappable_ratio'
		self.mappable_ratio_stat = 'mappable_ratio_graph'
		return self.mappable_ratio_stat

	def _redundant_ratio_info(self):
		""" Show redundant  ratio of the dataset in all historic data"""
		print 'redundant_ratio\n'
		self.redundant_ratio_stat = 'redundant_ratio_graph'
		return self.redundant_ratio_stat

	def _render(self):
		temp = self.template.render(MappingQC_check = self.MappingQC_check, basic_mapping_table = self.mappable_summary_stat, mappable_ratio_graph = self.mappable_ratio_stat, redundant_ratio_graph = self.redundant_ratio_stat)
		print temp
	def run(self):
		self.MappingQC_check = True
		""" Run some MappingQC function to get final result. """
		self.mappable_summary_stat = self._basic_mapping_statistics_info()
		self.mappable_ratio_stat = self._mappable_ratio_info()
		self.redundant_ratio_stat = self._redundant_ratio_info()
		self._render()
	def check():
		"""Check whether the MappingQC's result is ok. """
		print 'mapping qc pass or not'


class PeakcallingQC(QC_Controller):
	""" PeakcallingQC aims to describe the quality of peak calling result."""
	def __init__(self):
		super(PeakcallingQC, self).__init__()
		print 'init peak calling  qc'
	def _peak_summary_info(self):
		"""Basic statistic of peak calling result."""
		print ' summary_info'

	def _velcro_ratio_info(self):
		"""verlcro ratio is used to describe whether the peak is credible , The lower the result is more convenience.
		 The cumulative percentage plot can reflect the particularly dataset's verlcro ratio quality of all historic data."""
		print 'velcro_ratio_info '

	def _DHS_ratio_info(self):
		""" DHS ratio indicate the percentage of peaks overlap with DHSs site.
		The function can describe  the particularly dataset's DHS ratio quality of all historic data.
		"""
		print 'DHS_ratio_info'
	def _replicate_info(self):
		""" ReplicateQC aims to describe the similarity of replicate experiment. Venn diagram and correlation plot will be used."""
		print 'replicate_info\n'
	def _render(self):
		pass

	def run(self):
		""" Run some PeakcallingQC function to get final result. """
		self._peak_summary_info()
		self._velcro_ratio_info()
		self._DHS_ratio_info()
		self._replicate_info()
	def check():
		""" Check whether PeakcallingQC's result is ok. """
		print 'pass or not'

		
class AnnotationQC(QC_Controller):
	""" AnnotationQC aims to describe the quality of annotations after peak calling. """ 
	def __init__(self):
		super(AnnotationQC, self).__init__()
		print 'intialization of function qc'
	def _ceas_info(self):
		""" Describe peaks' distribution and relative position. """
		print 'ceas qc'
	def _conservation_info(self):
		""" Density plot of peaks's conservation."""
		print 'conservation qc'
	def _motif_info(self):
		""" QC of Sepose. """
		print 'motif info\n'
	def _render(self):
		pass
	def run(self):
		""" Run some AnnotationQC function. """
		self._ceas_info()
		self._conservation_info()
		self._motif_info()
	def check(self):
		""" Check whether AnnotationQC's result is ok. """
		print 'pass or not'

