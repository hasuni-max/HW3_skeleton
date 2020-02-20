import matplotlib.pyplot as plt

class roc(object):

	def __init__(self):
		self.tp = []
		self.fp = []

	def add_rates(self,truep,falsep):
		self.tp.append(truep)
		self.fp.append(falsep)

	def plot_ROC(self,lab):
		plt.plot(self.fp,self.tp,label=lab)
		plt.plot([0, 1], [0, 1], color='navy', linestyle='--')
		plt.ylabel('True Positive')
		plt.xlabel('False Positive')
		plt.title('Receiver operating characteristic')
		plt.legend(loc="lower right")
		#plt.show()	
	def new_curve(self):
		self.tp = []
		self.fp = []

	def show_plot(self):
		plt.show()

	def save_plot(self, filename):
		plt.savefig(filename)



if __name__ == "__main__":

	r = roc()

	r.plot_ROC()
	r.show_plot()

