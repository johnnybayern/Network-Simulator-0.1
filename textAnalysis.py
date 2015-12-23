'''
Created on May 15,2015
@author: Sagar Jha
'''

#from Tkinter import *
import nltk as nltk
from constants import *
from captureprint import capture_output


class AnalysisWidget(object):
    def __init__(self,master,text,textWidget,NB):
        top = self.top = Toplevel(master)

        self.textWidget = textWidget
        self.NB = NB

        self.text = text
        self.top.title("Text Analysis")
        # self.top.geometry('%dx%d+%d+%d' % (300,250,10,500))
        # self.top.resizable(False,False)
        self.top.resizable(TRUE,TRUE)

        self.inputFrame = LabelFrame(top, text="Keywords Search Functions")
        self.inputFrame.pack(padx=10, pady=10, expand=TRUE, fill=BOTH)

        self.plotFrame = LabelFrame(top, text="Plot Functions")
        self.plotFrame.pack(padx=10, pady=10, expand=TRUE, fill=BOTH)

        self.inputLabel = Label(self.inputFrame,text="Keywords")
        self.inputLabel.grid(row=0,column=0)
        # self.inputLabel.pack(side=LEFT)

        self.inputBox = Entry(self.inputFrame,font="Helvetica 12")
        self.inputBox.grid(row=0,column=1,columnspan=2)
        # self.inputbox.configure(fill=BOTH)
        # self.inputbox.pack(side=LEFT)
         
        self.searchText = Button(self.inputFrame,text="SearchText",width=25,command=self.__searchText)
        self.searchText.grid(row=1,column=1)
        # self.searchText.configure(fill=BOTH)
        # self.searchText.pack(side=BOTTOM)

        self.searchSText = Button(self.inputFrame,text="Search Similar Text",width=25,command=self.__searchSimilarText)
        self.searchSText.grid(row=2,column=1)
        # self.searchSText.config(fill=BOTH)
        # self.searchSText.pack(side=BOTTOM)

        self.searchCText = Button(self.inputFrame,text="Search Common Context Text", width=25, command=self.__seachCommonConText)
        self.searchCText.grid(row=3,column=1)
        # self.searchCText.configure(fill=BOTH)
        # self.searchCText.pack(fill=BOTH,expand=1)
        # self.searchCText.pack(side=BOTTOM)
        self.searchCollocations = Button(self.inputFrame, text="Search Collocations", width=25, command=self.__collocations)
        self.searchCollocations.grid(row=4,column=1)

        self.dPlot = Button(self.plotFrame,text="Dispersion Plot", width=25, command=self.__dispersionPlot)
        # self.dPlot.grid(row=1,column=0)
        # self.dPlot.configure(fill=BOTH)
        self.dPlot.pack(side=BOTTOM)
        
        self.fDist = Button(self.plotFrame,text="Frequency Distribution", width=25, command=self.__frequencyDistribution)
        # self.fDist.grid(row=2,column=0)
        # self.fDist.configure(fill=BOTH)
        self.fDist.pack(side=BOTTOM)


    def findValue(self):
        return self.Choice.get()
        
    def cleanup(self):
        self.top.destroy()

    def readTextbox(self):
        inputtext = self.inputBox.get()
        self.inputlist = []
       
        if inputtext != "":
            self.inputlist = inputtext.split(",")
        else:
            writeCalculations(self.textWidget,"Please enter words into the textbox!!!",True,self.NB)

    def __searchText(self):
        """

        :type self: object
        """
        # writeCalculations(self.textWidget,"-"*100 ,False,self.NB)
        print_seperator(self.textWidget,"Search Text" ,False,self.NB)
        # writeCalculations(self.textWidget,"-"*100 ,False,self.NB)
        self.readTextbox()
        for word in self.inputlist:
            values = capture_output(self.text.concordance, word)
            # values = values[:len(values)-1]
            # writeCalculations(self.textWidget,"*"*100 ,False,self.NB)
            # print_seperator(self.textWidget,"Displaying %d of %d matches for %s" % (len(values),len(values),word),False,self.NB)
            # writeCalculations(self.textWidget,"*"*100 ,False,self.NB)
            #for i in range(len(values)):
                # writeCalculations(self.textWidget,values[i][0]+" "+values[i][1]+" "+values[i][2],False,self.NB)
            writeCalculations(self.textWidget, values ,False,self.NB)

    def __searchSimilarText(self):
        # writeCalculations(self.textWidget,"-"*100 ,False,self.NB)
        print_seperator(self.textWidget,"Search Similar Text",False,self.NB)
        # writeCalculations(self.textWidget,"-"*100 ,False,self.NB)
        self.readTextbox()
        for word in self.inputlist:
            # writeCalculations(self.textWidget,"*"*100 ,False,self.NB)
            print_seperator(self.textWidget,"Similar text search for %s " % (word) ,False,self.NB)
            # writeCalculations(self.textWidget,"*"*100 ,False,self.NB)
            value = capture_output(self.text.similar, word)
            writeCalculations(self.textWidget,value,False,self.NB)

    def __seachCommonConText(self):
        # writeCalculations(self.textWidget,"-"*100 ,False,self.NB)
        print_seperator(self.textWidget,"Search common context",False,self.NB)
        # writeCalculations(self.textWidget,"-"*100 ,False,self.NB)
        self.readTextbox()
        value = capture_output(self.text.common_contexts, self.inputlist)
        writeCalculations(self.textWidget,value,False,self.NB)
       
    def __dispersionPlot(self):
        # writeCalculations(self.textWidget,"-"*100 ,False,self.NB)
        print_seperator(self.textWidget,"Dispersion Plot",False,self.NB)
        # writeCalculations(self.textWidget,"-"*100 ,False,self.NB)


        self.readTextbox()
        self.text.dispersion_plot(self.inputlist)

    def __frequencyDistribution(self):
        # writeCalculations(self.textWidget,"-"*100 ,False,self.NB)
        # writeCalculations(self.textWidget,"Frequency Distribution",False,self.NB)
        # writeCalculations(self.textWidget,"-"*100 ,False,self.NB)
        print_seperator(self.textWidget,"Frequency Distribution",False,self.NB)

        fdist1 = nltk.FreqDist(self.text)
        vocab1 = fdist1.keys()

        iNum = 50
        if len(vocab1) < iNum:
            iNum = len(vocab1)

        fdist1.plot(iNum,cumulative=True)

    def __collocations(self):
        print_seperator(self.textWidget,"Search Collocations",False,self.NB)
        value = capture_output(self.text.collocations)
        writeCalculations(self.textWidget,value,False,self.NB)