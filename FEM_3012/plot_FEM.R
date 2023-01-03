data=read.table("C:/Users/alexb/Desktop/FEM_project/FEM_3012/test.txt", sep=";")
x=as.numeric(gsub(",","",data[,2],fixed=TRUE))
print(x)
u=data[,1]
x=data[,2]
y=data[,3]
