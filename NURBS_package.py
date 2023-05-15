import math
import numpy as np

def B_spline_calculate_2d(knot_vector,order,vector):

    knot_vector = knot_vector_repeat_2d(knot_vector,order)

    len_knot = len(knot_vector)

    i = i_number(knot_vector,vector)

    B_spline = np.zeros(shape = (len_knot,order+1))

    B_spline[i][0] = 1

    for k in range(1,order+1):
            
        for n in range(i-k,i+1):

            if (knot_vector[n+k]-knot_vector[n] == 0 and (knot_vector[n+k+1]-knot_vector[n+1]) != 0):

                B_spline[n][k] = (knot_vector[n+k+1]-vector)/(knot_vector[n+k+1]-knot_vector[n+1])*B_spline[n+1][k-1]

            elif (knot_vector[n+k]-knot_vector[n] != 0 and (knot_vector[n+k+1]-knot_vector[n+1]) == 0):

                B_spline[n][k] = (vector-knot_vector[n])/(knot_vector[n+k]-knot_vector[n])*B_spline[n][k-1]

            elif (knot_vector[n+k]-knot_vector[n] == 0 and (knot_vector[n+k+1]-knot_vector[n+1]) == 0):

                B_spline[n][k] = 0

            else:

                B_spline[n][k] = (vector-knot_vector[n])/(knot_vector[n+k]-knot_vector[n])*B_spline[n][k-1]+(knot_vector[n+k+1]-vector)/(knot_vector[n+k+1]-knot_vector[n+1])*B_spline[n+1][k-1]

    return B_spline

def B_spline_point_calculate_2d(control_point,knot_vector,order,vector):

    size_control_point = control_point.shape
    len_control_point = size_control_point[1]

    i = i_number(knot_vector,vector)
    B_spline = B_spline_calculate_2d(knot_vector,order,vector)

    point = np.zeros(2)

    for d in range(2):

        point[d] = 0

        for m in range(order+1):

            if (i+m >= len_control_point):

                continue

            else:

                point[d] = point[d]+control_point[d][i+m]*B_spline[i+m][order]

    return point

def B_spline_point_derivative_calculate_2d(control_point,derivative_order,knot_vector,order,vector):

    size_control_point = control_point.shape
    len_control_point = size_control_point[1]

    control_point_derivative = control_point_derivative_calculate_2d(control_point,derivative_order,knot_vector,order)   
    B_spline = B_spline_calculate_2d(knot_vector,order,vector)

    B_spline_point_derivative = [0,0]

    for d in range(2):

        for j in range(len_control_point):

            B_spline_point_derivative[d] = B_spline_point_derivative[d]+control_point_derivative[d][j][derivative_order]*B_spline[j][order-derivative_order]

    return B_spline_point_derivative

def control_point_calculate_2d(data_point,first_tangent_vector,last_tangent_vector,knot_vector,order):

    size_data_point = data_point.shape
    len_data_point = size_data_point[1]

    control_point = np.zeros(shape = (2,len_data_point+2))

    max_knot_vector = max(knot_vector)    

    knot_vector = knot_vector_repeat_2d(knot_vector,order)

    for d in range(2):
        
        control_point[d][0] = data_point[d][0]
        control_point[d][1] = knot_vector[1+order]*first_tangent_vector[d]/order-control_point[d][0]
        control_point[d][len_data_point+1] = data_point[d][len_data_point-1]
        control_point[d][len_data_point] = control_point[d][len_data_point+1]-(max_knot_vector-knot_vector[len_data_point+order-2])*last_tangent_vector[d]/order

    data_point_matrix = np.delete(data_point,len_data_point-1,1)
    data_point_matrix = np.delete(data_point_matrix,0,1)

    size_data_point_matrix  = data_point_matrix.shape
    len_data_point_matrix  = size_data_point_matrix [1]

    control_point_matrix = np.zeros(shape = (2,len_data_point_matrix))
    coefficiency_matrix = np.zeros(shape = (len_data_point_matrix,len_data_point_matrix))

    for m in range(len_data_point_matrix):

        N = B_spline_calculate_2d(knot_vector,order,knot_vector[m+order+1])

        for d in range(2):

            if (m == 0):

                data_point_matrix[d][m] = data_point_matrix[d][m]-N[m+order+1][order]*control_point[d][1]

                for i in range(order-1):

                    coefficiency_matrix[m][m+i] = N[m+order+i+2][order]

            elif (m == len_data_point_matrix-1):

                data_point_matrix[d][m] = data_point_matrix[d][m]-N[m+order+3][order]*control_point[d][len_data_point]

                for i in range(order):

                    if (m+i-1 >= len_data_point_matrix):

                        continue

                    else:

                        coefficiency_matrix[m][m+i-1] = N[m+order+i+1][order]

            else:

                for i in range(order):

                    if (m+i-1 >= len_data_point_matrix):

                        continue

                    else:

                        coefficiency_matrix[m][m+i-1] = N[m+order+i+1][order]

    coefficiency_matrix_inv = np.linalg.inv(coefficiency_matrix)

    for d in range(2):
        
        control_point_matrix[d] = np.dot(coefficiency_matrix_inv,data_point_matrix[d])

        for m in range(2,len_data_point):

            control_point[d][m] = control_point_matrix[d][m-2]

    return control_point

def control_point_derivative_calculate_2d(control_point,derivative_order,knot_vector,order):

    size_control_point = control_point.shape
    len_control_point = size_control_point[1]

    control_point_derivative = np.zeros(shape = (2,len_control_point,derivative_order+1))

    knot_vector = knot_vector_repeat_2d(knot_vector,order)

    for d in range(2):

        for l in range(derivative_order+1):

            if (l == 0):

                for j in range(len_control_point-l):

                    control_point_derivative[d][j][l] = control_point[d][j]

            else:

                for j in range(len_control_point-l):

                    control_point_derivative[d][j][l] = (order-l+1)*(control_point_derivative[d][j+1][l-1]-control_point_derivative[d][j][l-1])/(knot_vector[j+order+1]-knot_vector[j+l])

    return control_point_derivative

def data_point_derivative_calculate_2d(data_point,derivative_order,knot_vector):

    size_data_point = data_point.shape
    len_data_point = size_data_point[1]

    data_point_derivative = np.zeros(shape = (2,len_data_point,derivative_order+1))

    knot_vector_derovative = knot_vector_derivative_calculate_2d(derivative_order,knot_vector)

    for d in range(2):

        for k in range(derivative_order+1):

            if (k == 0):

                for m in range(len_data_point-k):

                    data_point_derivative[d][m][k] = data_point[d][m]

            else:

                for m in range(len_data_point-k):

                    data_point_derivative[d][m][k] = (data_point_derivative[d][m+1][k-1]-data_point_derivative[d][m][k-1])/(knot_vector_derovative[m+1][k-1]-knot_vector_derovative[m][k-1])
    
    return data_point_derivative

def fit_2d(data_point,order,point_delta):

    size_data_point = data_point.shape
    len_data_point = size_data_point[1]

    vector_delta = point_delta/100

    knot_vector = initial_knot_vector_calculate_2d(data_point)

    max_knot_vector = max(knot_vector)

    control_point = initial_control_point_calculate_2d(data_point,knot_vector,order)

    B_spline_point_derivative = B_spline_point_derivative_calculate_2d(control_point,1,knot_vector,order,0)
    first_tangent_vector = [B_spline_point_derivative[0]/np.sqrt(B_spline_point_derivative[0]*B_spline_point_derivative[0]+B_spline_point_derivative[1]*B_spline_point_derivative[1]),B_spline_point_derivative[1]/np.sqrt(B_spline_point_derivative[0]*B_spline_point_derivative[0]+B_spline_point_derivative[1]*B_spline_point_derivative[1])]
    
    B_spline_point_derivative = B_spline_point_derivative_calculate_2d(control_point,1,knot_vector,order,max_knot_vector-vector_delta)
    last_tangent_vector = [B_spline_point_derivative[0]/np.sqrt(B_spline_point_derivative[0]*B_spline_point_derivative[0]+B_spline_point_derivative[1]*B_spline_point_derivative[1]),B_spline_point_derivative[1]/np.sqrt(B_spline_point_derivative[0]*B_spline_point_derivative[0]+B_spline_point_derivative[1]*B_spline_point_derivative[1])]

    knot_point_vector = knot_point_vector_calculate_2d(control_point,data_point,knot_vector,order,point_delta,vector_delta)

    knot_vector = knot_point_vector[0]
    point_vector = knot_point_vector[1]

    #point_vector = point_vector_calculate_2d(data_point,knot_vector,point_delta,vector_delta)

    max_knot_vector = max(knot_vector)

    control_point = control_point_calculate_2d(data_point,first_tangent_vector,last_tangent_vector,knot_vector,order)

    point_number = int((max(data_point[0])-min(data_point[0]))/point_delta)
    fit_point = np.zeros(shape = (2,point_number+1))

    for n in range(point_number+1):

        B_spline_point = B_spline_point_calculate_2d(control_point,knot_vector,order,point_vector[n])

        if (n == point_number):

            for d in range(2):

                fit_point[d][n] = data_point[d][len_data_point-1]

        else:

            for d in range(2):

                fit_point[d][n] = B_spline_point[d]

    return fit_point

def i_number(knot_vector,vector):

    len_knot = len(knot_vector)

    for n in range(len_knot-1):

        if ((vector-knot_vector[n]) >= 0 and (vector-knot_vector[n+1]) < 0):

            i = n
            break

    return i

def initial_control_point_calculate_2d(data_point,knot_vector,order):

    size_data_point = data_point.shape
    len_data_point = size_data_point[1]

    control_point = np.zeros(shape = (2,len_data_point+2))

    data_point_derivative = data_point_derivative_calculate_2d(data_point,1,knot_vector)

    first_tangent_vector = [data_point_derivative[0][0][1],data_point_derivative[1][0][1]]
    last_tangent_vector = [data_point_derivative[0][len_data_point-order+1][1],data_point_derivative[1][len_data_point-order+1][1]]

    max_knot_vector = max(knot_vector)

    knot_vector = knot_vector_repeat_2d(knot_vector,order)

    for d in range(2):
        
        control_point[d][0] = data_point[d][0]
        control_point[d][1] = knot_vector[1+order]*first_tangent_vector[d]/order-control_point[d][0]
        control_point[d][len_data_point+1] = data_point[d][len_data_point-1]
        control_point[d][len_data_point] = control_point[d][len_data_point+1]-(max_knot_vector-knot_vector[len_data_point+order-2])*last_tangent_vector[d]/order

    data_point_matrix = np.delete(data_point,len_data_point-1,1)
    data_point_matrix = np.delete(data_point_matrix,0,1)

    size_data_point_matrix  = data_point_matrix.shape
    len_data_point_matrix  = size_data_point_matrix [1]

    control_point_matrix = np.zeros(shape = (2,len_data_point_matrix))
    coefficiency_matrix = np.zeros(shape = (len_data_point_matrix,len_data_point_matrix))

    for m in range(len_data_point_matrix):

        N = B_spline_calculate_2d(knot_vector,order,knot_vector[m+order+1])

        for d in range(2):

            if (m == 0):

                data_point_matrix[d][m] = data_point_matrix[d][m]-N[m+order+1][order]*control_point[d][1]

                for i in range(order-1):

                    coefficiency_matrix[m][m+i] = N[m+order+i+2][order]

            elif (m == len_data_point_matrix-1):

                data_point_matrix[d][m] = data_point_matrix[d][m]-N[m+order+3][order]*control_point[d][len_data_point]

                for i in range(order):

                    if (m+i-1 >= len_data_point_matrix):

                        continue

                    else:

                        coefficiency_matrix[m][m+i-1] = N[m+order+i+1][order]

            else:

                for i in range(order):

                    if (m+i-1 >= len_data_point_matrix):

                        continue

                    else:

                        coefficiency_matrix[m][m+i-1] = N[m+order+i+1][order]

    coefficiency_matrix_inv = np.linalg.inv(coefficiency_matrix)

    for d in range(2):
        
        control_point_matrix[d] = np.dot(coefficiency_matrix_inv,data_point_matrix[d])

        for m in range(2,len_data_point):

            control_point[d][m] = control_point_matrix[d][m-2]

    return control_point

def initial_knot_vector_calculate_2d(data_point):

    size_data_point = data_point.shape
    len_data_point = size_data_point[1]

    knot_vector = np.zeros(len_data_point)

    for n in range(len_data_point-1):

        knot_vector[n+1] = knot_vector[n]+math.sqrt((data_point[0][n+1]-data_point[0][n])*(data_point[0][n+1]-data_point[0][n])+(data_point[1][n+1]-data_point[1][n])*(data_point[1][n+1]-data_point[1][n]))

    return knot_vector

def knot_point_vector_calculate_2d(control_point,data_point,knot_vector,order,point_delta,vector_delta):

    size_data_point = data_point.shape
    len_data_point = size_data_point[1]

    vector_number = int(max(knot_vector)/vector_delta)
    point_number = int((max(data_point[0])-min(data_point[0]))/point_delta)

    vector_length = 0
    CDF_point = 0
    m = 0

    new_knot_vector = np.zeros(len_data_point)
    point_vector = np.zeros(point_number+1)

    for u in range(vector_number+1):

        vector = u*vector_delta

        if (CDF_point >= m*point_delta):

            point_vector[m] = vector_length

            m = m+1

        B_spline_point_derivative = B_spline_point_derivative_calculate_2d(control_point,1,knot_vector,order,vector)

        i = i_number(knot_vector,vector)

        vector_length = vector_length+vector_delta*np.sqrt(1+(B_spline_point_derivative[1]*B_spline_point_derivative[1])/(B_spline_point_derivative[0]*B_spline_point_derivative[0]))/np.sqrt(1+(data_point[1][i+1]-data_point[1][i])*(data_point[1][i+1]-data_point[1][i])/((data_point[0][i+1]-data_point[0][i])*(data_point[0][i+1]-data_point[0][i])))
        CDF_point = CDF_point+vector_delta/np.sqrt(1+(data_point[1][i+1]-data_point[1][i])*(data_point[1][i+1]-data_point[1][i])/((data_point[0][i+1]-data_point[0][i])*(data_point[0][i+1]-data_point[0][i])))

        if (vector >= knot_vector[i]):

            new_knot_vector[i+1] = vector_length

    return new_knot_vector,point_vector

def knot_vector_calculate_2d_Hartley(control_point,data_point,order):

    size_control_point = control_point.shape
    len_control_point = size_control_point[1]

    control_point_length = np.zeros(len_control_point-1)

    for i in range(len_control_point-1):

        control_point_length[i] = np.sqrt((control_point[0][i+1]-control_point[0][i])*(control_point[0][i+1]-control_point[0][i])+(control_point[1][i+1]-control_point[1][i])*(control_point[1][i+1]-control_point[1][i]))

    size_data_point = data_point.shape
    len_data_point = size_data_point[1] 

    knot_vector  = np.zeros(len_data_point)
    knot_vector[len_data_point-1] = 1

    delta_length = np.zeros(len_data_point-1)
    total_length = 0

    for i in range(order,len_data_point+order-1):

        for j in range(i-order,i):

            if (j >= len_control_point-1):

                continue

            else:

                delta_length[i-order] = delta_length[i-order]+control_point_length[j]

                total_length = total_length+control_point_length[j]

    for i in range(1,len_data_point-1):
        
        knot_vector[i] = knot_vector[i-1]+delta_length[i-1]/total_length

    return knot_vector

def knot_vector_derivative_calculate_2d(derivative_order,knot_vector):

    len_knot_vector = len(knot_vector)

    knot_vector_derivative = np.zeros(shape = (len_knot_vector,derivative_order+1))

    for k in range(derivative_order+1):

        if (k == 0):
            
            for m in range(len_knot_vector-k):

                knot_vector_derivative[m][k] = knot_vector[m]

        else:

            for m in range(len_knot_vector-k):

                knot_vector_derivative[m][k] = (knot_vector_derivative[m][k-1]+knot_vector_derivative[m+1][k-1])/2

    return knot_vector_derivative

def knot_vector_repeat_2d(knot_vector,order):

    max_knot_vector = max(knot_vector)

    for k in range(order):
        
        knot_vector = np.insert(knot_vector,0,0)
        knot_vector = np.append(knot_vector,max_knot_vector)

    return knot_vector

def point_vector_calculate_2d(data_point,knot_vector,point_delta,vector_delta):

    vector_number = int(max(knot_vector)/vector_delta)
    point_number = int((max(data_point[0])-min(data_point[0]))/point_delta)

    vector_length = 0
    CDF_point = 0
    m = 0

    point_vector = np.zeros(point_number+1)

    for u in range(vector_number+1):

        vector = u*vector_delta

        if (CDF_point >= m*point_delta):

            point_vector[m] = vector_length

            m = m+1

        i = i_number(knot_vector,vector)

        vector_length = vector_length+vector_delta
        CDF_point = CDF_point+vector_delta/np.sqrt(1+(data_point[1][i+1]-data_point[1][i])*(data_point[1][i+1]-data_point[1][i])/((data_point[0][i+1]-data_point[0][i])*(data_point[0][i+1]-data_point[0][i])))

    return point_vector