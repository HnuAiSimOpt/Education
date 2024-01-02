f_inp = open("Job-2.inp", "r")
f_node = open("nodes.txt", "a+")
f_ele = open("elements.txt", "a+")
temp = f_inp.readline()
while temp != "*Node\n":
    temp = f_inp.readline()

while temp != "*Element, type=CPS3\n":
    temp = f_inp.readline()
    f_node.write(temp)

while temp != "*Nset, nset=Set-1, generate\n":
    temp = f_inp.readline()
    f_ele.write(temp)

f_node.close()
f_inp.close()
f_ele.close()
