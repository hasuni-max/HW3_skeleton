import matplotlib.pyplot as plt

class roc(object):

	"""
		Simple class for plotting mutliple curves on a single ROC graph. 

	"""

	def __init__(self):
		self.tp = []
		self.fp = []

	def add_rates(self,truep,falsep):
		"""
			Pass true positive and false positive rates to instance
		"""
		self.tp.append(truep)
		self.fp.append(falsep)

	def plot_ROC(self,lab):
		"""
			Once all tp and fp values have been appended to their respective lists, plot an ROC graph with a 
			given label, lab

		"""
		plt.plot(self.fp,self.tp,label=lab)
		plt.plot([0, 1], [0, 1], color='navy', linestyle='--')
		plt.ylabel('True Positive')
		plt.xlabel('False Positive')
		plt.title('Receiver operating characteristic')
		plt.legend(loc="lower right")
		#plt.show()	
	def new_curve(self):
		"""
			Once plot_ROC has been run, this method will empty tp and fp lists so that 
			a new curve can be added to the plot
		"""
		self.tp = []
		self.fp = []

	def show_plot(self):

		plt.show()

	def save_plot(self, filename):
		"""
			After all curves have been added to the ROC, save to filename
		"""
		plt.savefig(filename)



if __name__ == "__main__":

	r = roc()

	r.plot_ROC()
	r.show_plot()

