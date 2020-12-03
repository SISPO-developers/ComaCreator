import sys

if(len(sys.argv)>1):
    if(int(sys.argv[1])==0):
        import particle_sim 
    elif(int(sys.argv[1])==1):
        import convert_to_png
else:
    print("give me some params")
