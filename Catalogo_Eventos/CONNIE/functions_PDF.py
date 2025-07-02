import matplotlib.colors
import matplotlib.pyplot as plt
import matplotlib
import os

from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import A4, letter
from reportlab.lib.utils import ImageReader



def pdf_creator3x3(pdf_name, list_pixelmatrix):
    print('The PDF is being created')
    print(len(list_pixelmatrix), ' muons will be display')
    c = canvas.Canvas(pdf_name, letter)
    w, h = A4 # page size

    ### ===== Temporal Image ===== ###
    path_image = 'Perfil_Muon.jpg'
    color_matrix = [0.4, 0.4470, 0.2410]
    color_map ='jet'
    coord3x3 = [[0,0],[0,1], [0,2], [1,0], [1,1], [1,2], [2,0], [2,1], [2,2]] 


    n_last_events = 0
    n_muonstot = len(list_pixelmatrix)
    n_last_events = 0
    n_events = 0

    fig, axs = plt.subplots(ncols= 3, nrows= 3, figsize = [6,6], facecolor = color_matrix) ## Page will have 9 matrix 
    num_pages = 0


    if (n_muonstot%9) == 0:
        for index in range(0, len(list_pixelmatrix)):
            n_muonstot = n_muonstot - 1
            n_events = n_events + 1
            if n_events < 10: 
                # print('n_events: ', n_events)
                axs[coord3x3[n_events - 1][0],coord3x3[n_events - 1][1]].imshow(list_pixelmatrix[index], cmap = color_map)
                axs[coord3x3[n_events - 1][0],coord3x3[n_events - 1][1]].set_title('ID: ' + str(index))
                plt.tight_layout()

                if n_events == 9:
                    n_events = 10
            
            if n_events == 10:
                n_events = 0
                num_pages = num_pages + 1
                fig.savefig(path_image)
                img = ImageReader(path_image)

                img_w, img_h = img.getSize()    # width and high of the page

                c.drawImage(img, w - img_w + 10, h - img_h - 90)
                c.drawString(w/2, 40, str(num_pages))
                c.showPage()
                # plt.imshow(event)
                # plt.show()
                plt.close()
                os.system('rm Perfil_Muon.jpg') # if the provitional image name change then change this line

                fig, axs = plt.subplots(ncols= 3, nrows= 3, figsize = [6,6], facecolor = color_matrix)

    elif (n_muonstot%9) != 0:
        # print('No Hola')
        Resto = n_muonstot%9

        for index in range(0, len(list_pixelmatrix)):
            n_muonstot = n_muonstot - 1
            n_events = n_events + 1

            if n_muonstot >= Resto:
                if n_events < 10: 
                    # print('n_events: ', n_events)
                    axs[coord3x3[n_events - 1][0],coord3x3[n_events - 1][1]].imshow(list_pixelmatrix[index], cmap = color_map)
                    axs[coord3x3[n_events - 1][0],coord3x3[n_events - 1][1]].set_title('ID: ' + str(index))
                    plt.tight_layout()

                    if n_events == 9:
                        n_events = 10

                if n_events == 10:
                    n_events = 0
                    num_pages = num_pages + 1
                    fig.savefig(path_image)
                    img = ImageReader(path_image)

                    img_w, img_h = img.getSize()

                    c.drawImage(img, w - img_w+10, h - img_h-90)
                    c.drawString(w/2, 40, str(num_pages))
                    c.showPage()
                    # plt.imshow(event)
                    # plt.show()
                    plt.close()
                    os.system('rm Perfil_Muon.jpg')

                    fig, axs = plt.subplots(ncols= 3, nrows= 3, figsize = [6,6], facecolor = color_matrix)

            elif n_muonstot<Resto:
                n_last_events = n_last_events + 1

                axs[coord3x3[n_last_events - 1][0],coord3x3[n_last_events - 1][1]].imshow(list_pixelmatrix[index], cmap=color_map)
                axs[coord3x3[n_last_events - 1][0],coord3x3[n_last_events - 1][1]].set_title('ID: ' + str(index))
                plt.tight_layout()

                if n_muonstot == 0:
                    # plt.show()
                    num_pages = num_pages + 1
                    fig.savefig(path_image)
                    img = ImageReader(path_image)

                    img_w, img_h = img.getSize()

                    c.drawImage(img, w - img_w + 50, h - img_h - 60)
                    c.drawString(w/2, 40, str(num_pages))
                    c.showPage()
                    # plt.imshow(event)
                    # plt.show()
                    plt.close()
                    os.system('rm Perfil_Muon.jpg')

                    fig, axs = plt.subplots(ncols= 3, nrows= 3, figsize = [10,10], facecolor = color_matrix)
            
    plt.close()
    c.save()
    print('Done')

    return 0

def pdf_creator3x3_indexs(pdf_name, list_pixelmatrix, list_index):
    print('The PDF is being created')
    print(len(list_pixelmatrix), ' muons will be display')
    c = canvas.Canvas(pdf_name, letter)
    w, h = A4 # page size

    ### ===== Temporal Image ===== ###
    path_image = 'Perfil_Muon.jpg'
    color_matrix = [0.4, 0.4470, 0.2410]
    color_map ='jet'
    coord3x3 = [[0,0],[0,1], [0,2], [1,0], [1,1], [1,2], [2,0], [2,1], [2,2]] 
    Norm = matplotlib.colors.LogNorm(vmin=10**-1, vmax=10000 * 10**0)


    n_last_events = 0
    n_muonstot = len(list_pixelmatrix)
    n_last_events = 0
    n_events = 0

    fig, axs = plt.subplots(ncols= 3, nrows= 3, figsize = [6,6], facecolor = color_matrix) ## Page will have 9 matrix 
    num_pages = 0


    if (n_muonstot%9) == 0:
        for index in range(0, len(list_pixelmatrix)):
            n_muonstot = n_muonstot - 1
            n_events = n_events + 1
            if n_events < 10: 
                # print('n_events: ', n_events)
                axs[coord3x3[n_events - 1][0],coord3x3[n_events - 1][1]].imshow(list_pixelmatrix[index], cmap = color_map, norm = Norm)
                axs[coord3x3[n_events - 1][0],coord3x3[n_events - 1][1]].set_title('ID: ' + str(list_index[index]))
                plt.tight_layout()

                if n_events == 9:
                    n_events = 10
            
            if n_events == 10:
                n_events = 0
                num_pages = num_pages + 1
                fig.savefig(path_image)
                img = ImageReader(path_image)

                img_w, img_h = img.getSize()    # width and high of the page

                c.drawImage(img, w - img_w + 10, h - img_h - 90)
                c.drawString(w/2, 40, str(num_pages))
                c.showPage()
                # plt.imshow(event)
                # plt.show()
                plt.close()
                os.system('rm Perfil_Muon.jpg') # if the provitional image name change then change this line

                fig, axs = plt.subplots(ncols= 3, nrows= 3, figsize = [6,6], facecolor = color_matrix)

    elif (n_muonstot%9) != 0:
        # print('No Hola')
        Resto = n_muonstot%9

        for index in range(0, len(list_pixelmatrix)):
            n_muonstot = n_muonstot - 1
            n_events = n_events + 1

            if n_muonstot >= Resto:
                if n_events < 10: 
                    # print('n_events: ', n_events)
                    axs[coord3x3[n_events - 1][0],coord3x3[n_events - 1][1]].imshow(list_pixelmatrix[index], cmap = color_map, norm = Norm)
                    axs[coord3x3[n_events - 1][0],coord3x3[n_events - 1][1]].set_title('ID: ' + str(list_index[index]))
                    plt.tight_layout()

                    if n_events == 9:
                        n_events = 10

                if n_events == 10:
                    n_events = 0
                    num_pages = num_pages + 1
                    fig.savefig(path_image)
                    img = ImageReader(path_image)

                    img_w, img_h = img.getSize()

                    c.drawImage(img, w - img_w+10, h - img_h-90)
                    c.drawString(w/2, 40, str(num_pages))
                    c.showPage()
                    # plt.imshow(event)
                    # plt.show()
                    plt.close()
                    os.system('rm Perfil_Muon.jpg')

                    fig, axs = plt.subplots(ncols= 3, nrows= 3, figsize = [6,6], facecolor = color_matrix)

            elif n_muonstot<Resto:
                n_last_events = n_last_events + 1

                axs[coord3x3[n_last_events - 1][0],coord3x3[n_last_events - 1][1]].imshow(list_pixelmatrix[index], cmap=color_map, norm = Norm)
                axs[coord3x3[n_last_events - 1][0],coord3x3[n_last_events - 1][1]].set_title('ID: ' + str(list_index[index]))

                if n_muonstot == 0:
                    # plt.show()
                    num_pages = num_pages + 1
                    fig.savefig(path_image)
                    img = ImageReader(path_image)

                    img_w, img_h = img.getSize()

                    c.drawImage(img, w - img_w + 50, h - img_h - 60)
                    c.drawString(w/2, 40, str(num_pages))
                    c.showPage()
                    # plt.imshow(event)
                    # plt.show()
                    plt.close()
                    os.system('rm Perfil_Muon.jpg')

                    fig, axs = plt.subplots(ncols= 3, nrows= 3, figsize = [10,10], facecolor = color_matrix)
            
    plt.close()
    c.save()
    print('Done')

    return 0


def pdf_creator3x3_verticals(pdf_name, list_pixelmatrix):
    print('The PDF is being created')
    print(len(list_pixelmatrix), ' muons will be display')
    c = canvas.Canvas(pdf_name, letter)
    w, h = A4 # page size

    ### ===== Temporal Image ===== ###
    path_image = 'Perfil_Muon.jpg'
    color_matrix = [0.4, 0.4470, 0.2410]
    color_map ='jet'
    coord3x3 = [[0,0],[0,1], [0,2], [1,0], [1,1], [1,2], [2,0], [2,1], [2,2]] 

    Norm = matplotlib.colors.LogNorm(vmin=10**-1, vmax=5 * 10**0)


    n_last_events = 0
    n_muonstot = len(list_pixelmatrix)
    n_last_events = 0
    n_events = 0

    fig, axs = plt.subplots(ncols= 3, nrows= 3, figsize = [6,6], facecolor = color_matrix) ## Page will have 9 matrix 
    num_pages = 0


    if (n_muonstot%9) == 0:
        for index in range(0, len(list_pixelmatrix)):
            n_muonstot = n_muonstot - 1
            n_events = n_events + 1
            if n_events < 10: 
                # print('n_events: ', n_events)
                ylen, xlen = list_pixelmatrix[index].shape

                axs[coord3x3[n_events - 1][0],coord3x3[n_events - 1][1]].imshow(list_pixelmatrix[index], cmap = color_map, norm = Norm)
                axs[coord3x3[n_events - 1][0],coord3x3[n_events - 1][1]].set_title('ID: ' + str(index))
                axs[coord3x3[n_events - 1][0],coord3x3[n_events - 1][1]].set_ylabel('xlen: ' + xlen + ' , ylen: ' + ylen)
                plt.tight_layout()

                if n_events == 9:
                    n_events = 10
            
            if n_events == 10:
                n_events = 0
                num_pages = num_pages + 1
                fig.savefig(path_image)
                img = ImageReader(path_image)

                img_w, img_h = img.getSize()    # width and high of the page

                c.drawImage(img, w - img_w + 10, h - img_h - 90)
                c.drawString(w/2, 40, str(num_pages))
                c.showPage()
                # plt.imshow(event)
                # plt.show()
                plt.close()
                os.system('rm Perfil_Muon.jpg') # if the provitional image name change then change this line

                fig, axs = plt.subplots(ncols= 3, nrows= 3, figsize = [6,6], facecolor = color_matrix)

    elif (n_muonstot%9) != 0:
        # print('No Hola')
        Resto = n_muonstot%9

        for index in range(0, len(list_pixelmatrix)):
            n_muonstot = n_muonstot - 1
            n_events = n_events + 1

            if n_muonstot >= Resto:
                if n_events < 10: 
                    # print('n_events: ', n_events)
                    axs[coord3x3[n_events - 1][0],coord3x3[n_events - 1][1]].imshow(list_pixelmatrix[index], cmap = color_map, norm = Norm)
                    axs[coord3x3[n_events - 1][0],coord3x3[n_events - 1][1]].set_title('ID: ' + str(index))
                    plt.tight_layout()

                    if n_events == 9:
                        n_events = 10

                if n_events == 10:
                    n_events = 0
                    num_pages = num_pages + 1
                    fig.savefig(path_image)
                    img = ImageReader(path_image)

                    img_w, img_h = img.getSize()

                    c.drawImage(img, w - img_w+10, h - img_h-90)
                    c.drawString(w/2, 40, str(num_pages))
                    c.showPage()
                    # plt.imshow(event)
                    # plt.show()
                    plt.close()
                    os.system('rm Perfil_Muon.jpg')

                    fig, axs = plt.subplots(ncols= 3, nrows= 3, figsize = [6,6], facecolor = color_matrix)

            elif n_muonstot<Resto:
                n_last_events = n_last_events + 1

                axs[coord3x3[n_last_events - 1][0],coord3x3[n_last_events - 1][1]].imshow(list_pixelmatrix[index], cmap=color_map, norm = Norm)
                axs[coord3x3[n_last_events - 1][0],coord3x3[n_last_events - 1][1]].set_title('ID: ' + str(index))
                plt.tight_layout()

                if n_muonstot == 0:
                    # plt.show()
                    num_pages = num_pages + 1
                    fig.savefig(path_image)
                    img = ImageReader(path_image)

                    img_w, img_h = img.getSize()

                    c.drawImage(img, w - img_w + 50, h - img_h - 60)
                    c.drawString(w/2, 40, str(num_pages))
                    c.showPage()
                    # plt.imshow(event)
                    # plt.show()
                    plt.close()
                    os.system('rm Perfil_Muon.jpg')

                    fig, axs = plt.subplots(ncols= 3, nrows= 3, figsize = [10,10], facecolor = color_matrix)
            
    plt.close()
    c.save()
    print('Done')

    return 0

