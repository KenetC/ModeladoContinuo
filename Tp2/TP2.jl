
import Pkg 
Pkg.add("Images")
Pkg.add("FFTW")
Pkg.add("StatsBase")

using Images, FFTW, StatsBase
# completa la imagen de forma tal que sus dimensiones sean multiplos de 16
function completar_tam_16(im)
	tam1,tam2 = size(im)
	resto1 = tam1 % 16
	resto2 = tam2 % 16
	if (resto1 == 0 && resto2 == 0)
		return(im)
	else
		mitadizq = (16-resto2) ÷ 2
		mitadder= 16 - resto2 - mitadizq
		mitadarriba = (16-resto1) ÷ 2
		mitadabajo = 16 - resto1 - mitadarriba
		nueva_im = fill(RGB(0, 0, 0), tam1+16-resto1, tam2+16-resto2)
		nueva_im[mitadarriba+1:mitadarriba+tam1, mitadizq+1:mitadizq+tam2] = copy(im)
		return(nueva_im)
	end
end
# descompone la imagen a una matriz de tipo YCbCr luminosidad,Cb,Cr, -128 para centrala en 0 
function descomponer(im)
	im_16 = completar_tam_16(im)
	imagen1_ycbcr = YCbCr.(im_16)
	imagen1_channelview = channelview(imagen1_ycbcr)
	y = imagen1_channelview[1, :, :]
	cb = imagen1_channelview[2, :, :]
	nueva_im_cb = zeros(Float32, size(cb,1) ÷ 2, size(cb,2)÷2) 
	cr= imagen1_channelview[3, :, :]
	nueva_im_cr= zeros(Float32, size(cr,1) ÷ 2, size(cr,2)÷2) 
	
	for i in 1:(size(nueva_im_cb)[1])
		for j in 1:(size(nueva_im_cb)[2])
			nueva_im_cb[i,j]= (cb[2*i-1, 2*j-1] + cb[2*i, 2*j-1]+ cb[2*i-1, 2*j]+
									cb[2*i, 2*j])/4
			nueva_im_cr[i,j]= (cr[2*i-1, 2*j-1] + cr[2*i, 2*j-1]+ cr[2*i-1, 2*j]+
									cr[2*i, 2*j])/4
			
		end
	end
	return([y .- 128, nueva_im_cb .- 128, nueva_im_cr .- 128])
end
# Proceso inverso de la descomposicion 
function descomponer_inversa(im)
	y = im[1]
	cb = im[2]
	cr = im[3]
	nueva_im_cb = zeros(Float32, size(cb,1)*2, size(cb,2)*2)
	nueva_im_cr = zeros(Float32, size(cb,1)*2, size(cb,2)*2)
	
	for i in 1:size(cb,1)
	for j in 1:size(cb,2)
		nueva_im_cb[2*i-1,2*j-1]= cb[i, j]+ 128
		nueva_im_cb[2*i-1,2*j]= cb[i, j] + 128
		nueva_im_cb[2*i,2*j-1]= cb[i, j] + 128
		nueva_im_cb[2*i,2*j]= cb[i, j] +128
		
		nueva_im_cr[2*i-1,2*j-1]= cr[i, j]+128
		nueva_im_cr[2*i-1,2*j]= cr[i, j]+128
		nueva_im_cr[2*i,2*j-1]= cr[i, j]+128
		nueva_im_cr[2*i,2*j]= cr[i, j]+128
	end
	end
	nueva_im_y = y .+ 128
	nueva = zeros(Float32, 3, size(nueva_im_y,1), size(nueva_im_y,2))
	
	nueva[1,:, :] = nueva_im_y
	nueva[2,:, :] = nueva_im_cb
	nueva[3,:, :] = nueva_im_cr

	res = RGB.(colorview(YCbCr, nueva))
	return res
end
# Transformada 

function transformada(M)
	n,m = size(M)
	for i in 1:n÷8 
	for j in 1:m÷8
		vista = view(M, 8*(i-1)+1:8*i, 8*(j-1)+1:8*j)
		dct!(vista)
	end
	end
end
function inversa_transformada(M)
	n,m = size(M)
	for i in 1:n÷8 
	for j in 1:m÷8
		vista = view(M, 8*(i-1)+1:8*i, 8*(j-1)+1:8*j)
		idct!(vista)
	end
	end
end

# Cuantizacion
function cuantizar(M,quant::Matrix{UInt8})
	n,m = size(M)
	for i in 1:n÷8 
	for j in 1:m÷8
		vista0 = view(M, 8*(i-1)+1:8*i, 8*(j-1)+1:8*j)
		vista0 .= round.(vista0 ./ quant, digits=0)
	end
	end
	return convert(Matrix{Int8},M)
end
function inv_cuantizar(M,quant::Matrix{UInt8})
	n,m = size(M)
	for i in 1:n÷8 
	for j in 1:m÷8
		vista = view(M, 8*(i-1)+1:8*i, 8*(j-1)+1:8*j)
		vista .= vista .* quant
	end
	end
end

# Matrices de cuantizacion que vamos a usar
# la primera matriz conserva mejor la calidad de la imagen
# la segunda baja mas la calidad, pues hace que perdamos precision en los valores en la compresion
# esto es producto del redondeo, cuando dividimos casillero a casillero los valores de la matriz de la imagen, con estas matrices.
function quant1()
	return convert(Matrix{UInt8},[16 11 10 16 24 40 51 61;
	12 12 14 19 26 58 60 55;
	14 13 16 24 40 57 69 56;
	14 17 22 29 51 87 80 62;
	18 22 37 56 68 109 103 77;
	24 35 55 64 81 104 113 92;
	49 64 78 87 103 121 120 101;
	72 92 95 98 112 100 103 99])
end
function quant2()
	return convert(Matrix{UInt8},[ 80  60  50  80  120 200 255 255;
	55  60  70  95  130 255 255 255;
	70  65  80  120 200 255 255 255;
	70  85  110 145 255 255 255 255;
	90  110 185 255 255 255 255 255;
	120 175 255 255 255 255 255 255;
	245 255 255 255 255 255 255 255;
	255 255 255 255 255 255 255 255])
end

#Compresion 
# Para el mejor manejo del orden zigzag, generamos los indices
function orden_zigzag()
	res = fill([UInt8(0),UInt8(0)],64)
	n = 8
	iterator = 1 
	for k in 2:(n + n - 1)
		if k % 2 == 0
			for i in max(1, k - n):min(k - 1, n)
				res[iterator] = [UInt8(k - i),UInt8(i)]
				iterator = iterator + 1 
			end
		else
			for i in max(1, k - n):min(k - 1, n)
				res[iterator] = [UInt8(i),UInt8(k-i)]
				iterator = iterator + 1
			end
		end
	end
	res[iterator] = [UInt8(8),UInt8(8)]
	return res
end

function compresion_matriz(M)
	n,m = size(M)
	res = zeros(Int8,n*m)
	orden_zigzag1 = orden_zigzag()
	ind = 1
	size_total = 1
	for i in 1:n÷8
	for j in 1:m÷8
		vista = view(M, 8*(i-1)+1:8*i, 8*(j-1)+1:8*j)
		vect = zeros(Int8,64)
		for i in 1:64
			vect[i] = vista[orden_zigzag1[i][1],orden_zigzag1[i][2]]
		end
	
		vals, reps = rle(vect)
		v = vcat(reps,vals)
		
		res[ind:size(v,1)+ind-1] = v
		
		size_total = size_total + size(v,1)
		ind = ind + size(v,1)  
	end
	end
	return res[1:size_total-1]
end 

function inversa_compresion(A,bloques_fils,bloques_cols)
	M = zeros(bloques_fils*8,bloques_cols*8)
	orden_zigzag1 = orden_zigzag()
	i = 1
	j = 1
	for ind in 1:size(A,1)
		if ind % bloques_cols == 0 
			vista = view(M, 8*(i-1)+1:8*i, 8*(j-1)+1:8*j)
			V = inverse_rle(A[ind][1],A[ind][2])
			for i in 1:64 
				vista[orden_zigzag1[i][1],orden_zigzag1[i][2]] = V[i]
			end
			i += 1
			j = 1
			#println((i,j))
		else 
			vista = view(M, 8*(i-1)+1:8*i, 8*(j-1)+1:8*j)
			V = inverse_rle(A[ind][1],A[ind][2])	
			for i in 1:64 
				vista[orden_zigzag1[i][1],orden_zigzag1[i][2]] = V[i]
			end
			
			#println(vista)
			j += 1
			#println((i,j))
		end
	end
	return M
end


# Lectura y guardado 
# recibe, io que es el objeto que usamos para escribir o leer un archivo, en este caso
# vamos a leer una matriz, para ellos vamos a necesitar las dimensiones del mismo.
function leer_matriz(bloq_cols,bloq_fils,io)
	i = 1 
	j = 1 
	termino = false
	suma_parcial = 0 
	s = 0 
	reps = zeros(Int8,0)
	vals = zeros(Int8,0)
	A = []
	while(!termino)
		if suma_parcial == 64 
			while s != 0 
				val = read(io,Int8)
				push!(vals,val)
				s -= 1
			end
			suma_parcial = 0 
			
			push!(A,[vals,reps])
			reps = zeros(Int8,0)
			vals = zeros(Int8,0)
			if(j == bloq_cols && i < bloq_fils)
				j = 1
				i += 1
			elseif(j == bloq_cols && i == bloq_fils)
				termino = true
			else
				j += 1 
			end 
		else 
			rep = read(io,Int8)
			
			suma_parcial += rep 
			push!(reps,rep)
			s += 1
		end
	end
	return A 
end
# aca leemos todo el archivo para luego proceder con la descompresion.
function leer_EXT(nombre::String)
	io = open(nombre*".ext", "r")
	n = read(io, UInt16)
	m = read(io, UInt16)
	quant = zeros(UInt8,8,8)
	for i in 1:8
	for j in 1:8
		quant[i,j] = read(io,UInt8)
	end
	end
	
	bloq_fils = n ÷ 8
	bloq_cols = m ÷ 8
	bloq_cols2 = bloq_cols ÷ 2
	bloq_fils2 = bloq_fils ÷ 2
	
	A = leer_matriz(bloq_cols,bloq_fils,io)
	M1 = inversa_compresion(A,bloq_fils,bloq_cols)
	B = leer_matriz(bloq_cols2,bloq_fils2,io)		
	M2 = inversa_compresion(B,bloq_fils2,bloq_cols2)
	C = leer_matriz(bloq_cols2,bloq_fils2,io)		
	M3 = inversa_compresion(C,bloq_fils2,bloq_cols2)
	close(io)
	return quant,[M1,M2,M3]
end
# guardamos en un archivo la compresion hecha.
function guardar_en_EXT(M,n,m,nombre,quant::Matrix{UInt8})
	io = open(nombre*".ext","w")
	write(io,n)
	write(io,m)
	for i in 1:8
	for j in 1:8
		write(io,quant[i,j])
	end
	end
	for a in M[1]
		write(io,a)
	end		
	for a in M[2]
		write(io,a)
	end
	for a in M[3]
		write(io,a)
	end
	close(io)
end
# Parte final 

# Recibe la ruta a la imagen, matriz de cuantizacion a usar para la compresion.
function primer_paso(path_img,quant::Matrix{UInt8})
	img = load(path_img)
	img_descomp = descomponer(img)
	#println("Descomponer: ")
	#println(img_descomp[1])
	#println(img_descomp[2])
	#println(img_descomp[3])
	
	#Transformamos
	transformada(img_descomp[1])
	transformada(img_descomp[2])
	transformada(img_descomp[3])

	#println("transformada: ")
	#println(img_descomp[1])
	#println(img_descomp[2])
	
	img_cuantizada1 = cuantizar(img_descomp[1],quant)
	img_cuantizada2 = cuantizar(img_descomp[2],quant)
	img_cuantizada3 = cuantizar(img_descomp[3],quant)
	#println("Cuantizar: ")
	#println(img_cuantizada1)
	n,m = size(img_cuantizada1)

	y = compresion_matriz(img_cuantizada1)
	cb = compresion_matriz(img_cuantizada2)
	cr = compresion_matriz(img_cuantizada3)
	#println("Compresion: ")
	#println(y)
	
	nombre = path_img[begin:end-4]
	guardar_en_EXT([y,cb,cr],convert(UInt16,n),convert(UInt16,m),nombre,quant)
	return nombre
end
function segundo_paso(nombre)
	# nos devuelve array de 3 matrices, intensidad, cb y cr 
	quant,M = leer_EXT(nombre)
	# inversa de cuantizacion 
	inv_cuantizar(M[1],quant)
	inv_cuantizar(M[2],quant)
	inv_cuantizar(M[3],quant)

	#INVERSA DE TRANSFORMADA 	
	inversa_transformada(M[1])
	inversa_transformada(M[2])
	inversa_transformada(M[3])

	img = descomponer_inversa(M)
	#println(typeof(img))
	return img
end


#Probamos con "bolitas.bmp" y "paisaje.bmp" y anduvo relativamente rapido para cualquier matriz de cuantizacion 
#, sin embargo para "Meisje_met_de_parel.jpg" tardaba unos segundos mas, que es lo que se espera por el tamanio de la matriz de la imagen, para ello usamos la matriz de quant3 que es la que 
# hace perder mas calidad en la imagen, sin embargo no hizo un cambio significativo en terminos del tiempo.

# la ruta era simple pues tenia la imagen en la misma ruta del codigo .jl por ejemplo
# nombre1 = primer_paso("paisaje.bmp",quant1()) o podemos usar quant2

# Uso el nombre para ubicar el archivo comprimido para hacer el proceso inverso
# segundo_paso(nombre1)

# Para la imagen mas "grande"

# nombre2 = primer_paso("Meisje_met_de_parel.jpg",quant2())

# segundo_paso(nombre2)

