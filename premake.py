import os

def get_fotran_files():
  with open('files','w') as outfile:
    outfile.write("files=")
    for root, subdir, allfiles in os.walk("src/"):
      for each in allfiles:
        string=os.path.join(root,each)
        if string.endswith(".f90"):
          print>>outfile,string[4:]+" \\"

def create_subfolder_in_obj():
    for root, subdir, allfiles in os.walk("src/"):
      # creating directory tree in object folder 
      # for object being created by make
      for each in root:
        mypath = "obj/"+root[4:]
        if not os.path.isdir(mypath):
          print mypath
          os.makedirs(mypath)


get_fotran_files()
create_subfolder_in_obj()
