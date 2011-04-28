#Evaluates the sample JSON file.
#Use this with python versions < 2.6
#JSON files have a similar structure to python Dict lists
#So it is possible to just eval the contents as if it is python code.
#Not pretty but gets the job done

def ReadJSON(fn):
  f = open(fn,"r")
  input = f.read()
  js = eval(input)
  return js
