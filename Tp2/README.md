Tenemos 2 versiones, Tp2.jl es un notebook, donde podemos interactuar de forma mas comoda con el codigo, TP2.jl es un script que podemos correr desde la terminal para directamente usar la compresion.

Pasos a seguir para usar el codigo desde la terminal: 

* Abrimos la terminal

* Corremos Julia: Julia 

* include("ruta/TP2.jl")

Podemos usar diferentes matrices de cuantizacion tenemos quant1() y quant2():

Compresion
* nombre = primer_paso("bolitas.bmp",quant1()) 

Descompresion
* imagen = segundo_paso(nombre)

