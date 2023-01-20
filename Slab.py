# this program allows to size the slab panels according to the BAEL 91 modify 99

# By KOUDAMA HOSE WISDOM Civil engineering student
# mail : wisdomkoudama@gmail.com

import streamlit as st
import math
import pandas as pd
import numpy as np
import cv2

st.markdown("by wisdomkoudama@gmail.com")
data = {'phi': [1,2,3,4,5,6,7,8,9,10],
        '5': [0.2,0.39,0.59,0.79,0.98,1.18,1.37,1.57,1.77,1.96],
        '6': [0.28,0.57,0.85,1.13,1.41,1.7,1.98,2.26,2.54,2.83],
        '8': [0.5,0.57,1.51,2.01,2.51,3.02,3.52,4.02,4.52,5.03],
        '10': [0.79,1.57,2.36,3.14,3.93,4.71,5.5,6.28,7.07,7.85],
        '12': [1.13,2.26,3.39,4.52,5.65,6.79,7.92,9.05,10.18,11.31],
        '14': [1.54,3.08,4.62,6.16,7.7,9.24,10.78,12.32,13.85,15.39],
        '16': [2.01,4.02,6.03,8.04,10.05,12.06,14.07,16.08,18.1,20.11],
        '20': [3.14,6.28,9.42,12.57,15.71,18.85,21.99,15.13,28.27,31.42],
        '25': [4.91,9.82,14.73,14.73,19.54,29.45,34.36,39.27,44.18,49.09],
        '32': [8.04,16.08,24.13,32.17,40.21,48.25,56.3,64.34,72.38,80.42,],
        '40': [12.57,25.13,37.7,50.27,62.83,75.4,87.96,100.5,113.1,125.7],
        }
df = pd.DataFrame(data = data)

st.markdown("""
<style>
#MainMenu{
    visibility : hidden;
}
</>
""", unsafe_allow_html= True)

st.title("SLAB PANEL")
st.markdown("---")
st.markdown("### this program allows to size the solid slab panels according to the BAEL 91 Modified 99")
col1, col2 = st.columns(2)
projectName = col1.text_input("Project Name")
informations = col2.text_input("info")
st.dataframe(df)

with st.sidebar:
    st.markdown("# parameters")
    LX = float(st.text_input(" Lx small side [m]", 1))
    LY = float(st.text_input(" Ly large side [m]", 1))
    G = float(st.text_input("Permanent loads G [Kn/m²]", 1))
    Q = float(st.text_input("exploitation load Q [Kn/m²]", 1))
    col1, col2 = st.columns(2)
    acier = col1.radio("Type of   steel", options=("FE400", "FE500"))
    fisuration =col2.radio("Cracking ", options= ("FPP", "FP/FTP"))
    fc = float(st.selectbox("Fc28 [MPa]", options= (16, 20, 22, 25, 30)))



#alpha pout voir s'il s'agit d'une dalle ou une poutre dalle
    alpha = round((LX / LY),2)
#------------------------------------combinaison de charge----------------------------------------------------------------
    pu = 1.35 * G + 1.5 * Q
    pser = G + Q

#---------------------------------fonction de calcul du moment pour alpha < 0,4--------------------------------------------
    def moment(x):
        return x*math.pow(LX, 2)/8
# --------------------------------------définition de Ux et Uy ---------------------------------------------------
    Ux = 1/(1 + 2.4 * math.pow(alpha, 3))
    Uy = math.pow(alpha,3) * (1.9 - 0.9 * alpha)
#---------------------------------------fonction-------------------------------------------------------------
    def sectionELU(fc, Mu, b, h):
        try:
            fbu = (0.85 * fc) / 1.5
            NUu = (Mu) / (b * math.pow((0.9 * h), 2) * fbu)
            if NUu <= 0.391:
                alphaU = 1.25 * (1 - math.sqrt(1 - (2 * NUu)))
                st = (3.5 * (1 - alphaU)) / alphaU
                if st <= 1.739:
                    contrainte = 200000 * st
                    Asu = Mu / (contrainte * (0.9 * h) * (1 - 0.4 * alphaU))
                else:
                    contrainte = 400 / 1.15
                    Asu = Mu / (contrainte * (0.9 * h) * (1 - (0.4 * alphaU)))
                return float(Asu)
            else:
                alphaU = Mu - (0.391 * b * math.pow((0.9 * h), 2) * fbu)
                if alphaU >= 0.4 * Mu:
                    MUcalcul = 0.6 * NUu
                else:
                    MUcalcul = 0.391
                alphacalcul = 1.25 * (1 - math.sqrt(1 - (2 * MUcalcul)))
                epsiprimSC = 3.5 * (1 - ((0.1 * h) / (alphacalcul * (0.9 * h))))
                epsiST = 3.5 * ((1 - alphacalcul) / alphacalcul)
                if epsiprimSC <= 1.739:
                    contraintreprimSC = 200000 * epsiprimSC
                else:
                    contraintreprimSC = 400 / 1.15
                if epsiST <= 1.739:
                    contraintST = 200000 * epsiST
                else:
                    contraintST = 400 / 1.15
                AprimSU = (Mu - MUcalcul * b * math.pow((0.9 * h), 2) * fbu) / (
                        contraintreprimSC * ((0.9 * h) - (0.1 * h)))
                ASU = ((AprimSU * contraintreprimSC) + (
                        0.8 * alphacalcul * b * (0.9 * h) * fbu)) / contraintST
                return float(ASU)
        except:
            pass
#--------------------------------------fonction for FP or FTP------------------------------------------------
    def sectionELS(fc1, Mser, b1, h1):
        try:
            alphaLim = (9 * fc1) / ((9 * fc1) + 200)
            MlimSER = 0.1 * fc1 * b1 * math.pow((0.9 * h1), 2) * alphaLim * (3 - alphaLim)
            if Mser <= MlimSER:
                Lamda = 1 + (30 * Mser) / (b1 * math.pow((0.9 * h1), 2) * 200)
                mediaire = math.pow(Lamda, (-3 / 2))
                phy = math.degrees(math.acos(mediaire))
                alpha = 1 + (2 * math.sqrt(Lamda) * math.cos(240 + phy / 3))
                Asser = Mser / (200 * (0.9 * h1) * (1 - (alpha / 3)))
            else:
                Lamda = Mser - MlimSER
                contraintPrimSC = (9 * fc1) * (1 - (0.1 * h1) / (alphaLim * (0.9 * h1)))
                AprimSser = Lamda / (contraintPrimSC * ((0.9 * h1) - (0.1 * h1)))
                Asser = ((AprimSser * contraintPrimSC) + (0.3 * alphaLim * b1 * (0.9 * h1) * fc1)) / 200
            return float(Asser)
        except:
            pass

#--------------------------------------fonction pour acier min--------------------------------------------------
    def secton_retenu(x, y):
        if x < y:
            x = y
        return x
#--------------------------------------détermination selon type de dalle-------------------------------------
    if alpha < 0.4:
        hauteurdalle = LX/20.0
        section_min_y = 8 * hauteurdalle
        section_min_x = section_min_y * (3 - alpha) / 2
        hauteurdalle *= 1000  # Conversion [eng]
        if fisuration == "FPP":
            Mox = moment(pu)
            Mox *= 1000000 # Conversion [eng]
            Moy = 0
            section_acier = float(sectionELU(fc, Mox, 1000, hauteurdalle))
            section_acier = secton_retenu(section_acier, section_min_y)
            section_acier_x = section_acier / 400
            espacement_x = min((3 * hauteurdalle), 33)
            espacement_y = min((4 * hauteurdalle), 45)
        elif fisuration == "FP/FTP":
            Mox = moment(pser)
            Mox *= 1000000  # Conversion [eng]
            section_acier = sectionELS(fc, Mox, 1000, hauteurdalle)
            section_acier = secton_retenu(section_acier, section_min_y)
            section_acier_x = section_acier / 400
            espacement_x = min((1.5 * hauteurdalle), 20)
            espacement_y = min((1.5 * hauteurdalle), 20)
    else:
        hauteurdalle = LX/30.0
        section_min_y = 8 * hauteurdalle
        section_min_x = section_min_y * (3 - alpha) / 2
        Mox = Ux * pu * math.pow(LX, 2)
        Moy = Uy * Mox

        if fisuration == "FPP":
            section_acier_x = sectionELU(fc, Mox, 1000, hauteurdalle)
            section_acier_x = secton_retenu(section_acier_x, section_min_x)
            section_acier_y = sectionELU(fc, Moy,1000, hauteurdalle)
            section_acier_y = secton_retenu(section_acier_y, section_min_y)
            espacement_x = min((3 * hauteurdalle), 33)
            espacement_y = min((4 * hauteurdalle), 45)
        elif fisuration == "FP/FTP":
            section_acier_x = sectionELS(fc, Mox, 1000, hauteurdalle)
            section_acier_x = section_acier_x/100
            section_acier_x  = secton_retenu(section_acier_x, section_min_x)
            section_acier_y = sectionELS(fc, Moy, 1000, hauteurdalle)
            section_acier_y = section_acier_y/100
            section_acier_y = secton_retenu(section_acier_y, section_min_y)
            espacement_x = min((1.5 * hauteurdalle), 20)
            espacement_y = min((1.5 * hauteurdalle), 20)


    def resultat_acier():
        if alpha < 0.4:
            col5, col6, col7= st.columns(3)
            col5.write("steel section along the y-axis")
            col6.write(round(section_acier/100, 2))
            col7.write(" Cm²")
            col5.write("Steel section along the x-axis ")
            section_x_retenu = secton_retenu(section_acier/400, section_min_x)
            col6.write(round(section_x_retenu, 2))
            col7.write(" Cm²")
        else:
            print(type(section_acier_x))
            col8, col9, col10 = st.columns(3)
            col8.write("Steel section along the y-axis ")
            col9.write(round(section_acier_y,2) )
            col10.write(" Cm²")
            col8.write("Steel section along the x-axis ")
            col9.write(round(section_acier_x, 2) )
            col10.write(" Cm²")


    resultat_0 = st.sidebar.button("steel section", on_click = resultat_acier)
    st.markdown(" # Choice of steel")
    col3, col4 = st.columns(2)
    nombre_y = int(col3.text_input("N-Y", 1))
    acier_y = int(col4.text_input("HA-Y", 1))
    nombre_x = int(col3.text_input("N-X", 1))
    acier_x = int(col4.text_input("HA-X", 1))

    img = np.ones((1000, 1000, 3), dtype=np.uint8)
    img = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
    cv2.putText(img, str(nombre_y), (50, 385), cv2.FONT_HERSHEY_COMPLEX, 0.7, (255, 255, 255), 1)
    cv2.putText(img, "HA", (70, 385), cv2.FONT_HERSHEY_COMPLEX, 0.7, (255, 255, 255), 1)
    cv2.putText(img, str(acier_x), (108, 385), cv2.FONT_HERSHEY_COMPLEX, 0.7, (255, 255, 255), 1)

    cv2.putText(img, str(nombre_x), (10, 430), cv2.FONT_HERSHEY_COMPLEX, 0.7, (255, 255, 255), 1)
    cv2.putText(img, "HA", (10, 460), cv2.FONT_HERSHEY_COMPLEX, 0.7, (255, 255, 255), 1)
    cv2.putText(img, str(acier_x), (10, 490), cv2.FONT_HERSHEY_COMPLEX, 0.7, (255, 255, 255), 1)
    point = (int(LY) * 100) + 50
    cv2.rectangle(img, (50, 280), (point, 350), (255, 255, 255), 2)
    cv2.line(img, (60, 300), (((int(LY) * 100) + 30), 300), (123, 0, 255), 2)
    for i in range(nombre_y):
        cv2.circle(img, ((80 + (i * (int(LY) * 100 // nombre_y))), 310), 5, (255, 255, 255), -1)
        cv2.circle(img, ((80 + (i * (int(LY )* 100 // nombre_y))), 340), 5, (255, 255, 255), -1)
    cv2.line(img, (60, 330), (((int(LY)* 100) + 30), 330), (123, 0, 255), 2)
    cv2.rectangle(img, (50, 400), ((50 + int(LY) * 100), (400 + (int(LX) * 100))), (255, 255, 255), 2)
    for i in range(nombre_x):
        cv2.line(img, (60, (430 + (i * 10))), ((int(LY)* 100), (430 + (i * 10))), (123, 0, 255), 2)  # h
    for i in range(nombre_y):
        cv2.line(img, ((70 + (i * 10)), 410), ((70 + (i * 10)), ((int(LX)*100) + 380)), (123, 0, 255), 2)
    cv2.putText(img, "Calculation note", (290, 30), cv2.FONT_HERSHEY_COMPLEX, 1, (255, 255, 255), 1)
    cv2.line(img, (290, 35), (580, 35), (255, 255, 255), 2)

    cv2.putText(img, "Project Name :", (80, 65), cv2.FONT_HERSHEY_COMPLEX, 0.7, (255, 255, 255), 1)
    cv2.putText(img, str(projectName), (400, 65), cv2.FONT_HERSHEY_COMPLEX, 0.7, (255, 255, 255), 1)
    cv2.putText(img, "Info :", (196, 100), cv2.FONT_HERSHEY_COMPLEX, 0.7, (255, 255, 255), 1)
    cv2.putText(img, str(informations), (400, 100), cv2.FONT_HERSHEY_COMPLEX, 0.7, (255, 255, 255), 1)

    cv2.putText(img, "Lx:", (50, 130), cv2.FONT_HERSHEY_COMPLEX, 0.7, (255, 255, 255), 1)
    cv2.putText(img, str(LX), (100, 130), cv2.FONT_HERSHEY_COMPLEX, 0.7, (255, 255, 255), 1)
    cv2.putText(img, "Ly :", (200, 130), cv2.FONT_HERSHEY_COMPLEX, 0.7, (255, 255, 255), 1)
    cv2.putText(img, str(LY), (250, 130), cv2.FONT_HERSHEY_COMPLEX, 0.7, (255, 255, 255), 1)
    cv2.putText(img, "G :", (400, 130), cv2.FONT_HERSHEY_COMPLEX, 0.7, (255, 255, 255), 1)
    cv2.putText(img, str(G), (450, 130), cv2.FONT_HERSHEY_COMPLEX, 0.7, (255, 255, 255), 1)
    cv2.putText(img, "Q :", (550, 130), cv2.FONT_HERSHEY_COMPLEX, 0.7, (255, 255, 255), 1)
    cv2.putText(img, str(Q), (600, 130), cv2.FONT_HERSHEY_COMPLEX, 0.7, (255, 255, 255), 1)

    cv2.putText(img, "alpha:", (50, 160), cv2.FONT_HERSHEY_COMPLEX, 0.7, (255, 255, 255), 1)
    cv2.putText(img, str(alpha), (140, 160), cv2.FONT_HERSHEY_COMPLEX, 0.7, (255, 255, 255), 1)
    cv2.putText(img, "Pu :", (240, 160), cv2.FONT_HERSHEY_COMPLEX, 0.7, (255, 255, 255), 1)
    cv2.putText(img, str(round(pu, 2)), (290, 160), cv2.FONT_HERSHEY_COMPLEX, 0.7, (255, 255, 255), 1)
    cv2.putText(img, "Pser :", (400, 160), cv2.FONT_HERSHEY_COMPLEX, 0.7, (255, 255, 255), 1)
    cv2.putText(img, str(round(pser,2)), (500, 160), cv2.FONT_HERSHEY_COMPLEX, 0.7, (255, 255, 255), 1)

    if alpha < 0.4:
        cv2.putText(img, "Moment:", (10, 190), cv2.FONT_HERSHEY_COMPLEX, 0.7, (255, 255, 255), 1)
        cv2.putText(img, str(round(Mox, 2)), (110, 190), cv2.FONT_HERSHEY_COMPLEX, 0.7, (255, 255, 255), 1)
        cv2.putText(img, "x-steel :", (300, 190), cv2.FONT_HERSHEY_COMPLEX, 0.7, (255, 255, 255), 1)
        cv2.putText(img, str(round(section_acier_x, 2)), (420, 190), cv2.FONT_HERSHEY_COMPLEX, 0.7, (255, 255, 255), 1)
        cv2.putText(img, " y-steel :", (600, 190), cv2.FONT_HERSHEY_COMPLEX, 0.7, (255, 255, 255), 1)
        cv2.putText(img, str(round((section_acier/100), 2)), (730, 190), cv2.FONT_HERSHEY_COMPLEX, 0.7, (255, 255, 255), 1)
    else:
        cv2.putText(img, "Moment_x :", (50, 190), cv2.FONT_HERSHEY_COMPLEX, 0.7, (255, 255, 255), 1)
        cv2.putText(img, str(round(Mox, 2)), (200, 190), cv2.FONT_HERSHEY_COMPLEX, 0.7, (255, 255, 255), 1)
        cv2.putText(img, "Moment_y :", (360, 190), cv2.FONT_HERSHEY_COMPLEX, 0.7, (255, 255, 255), 1)
        cv2.putText(img, str(round(Moy, 2)), (520, 190), cv2.FONT_HERSHEY_COMPLEX, 0.7, (255, 255, 255), 1)
        cv2.putText(img, "x-steel :", (50, 220), cv2.FONT_HERSHEY_COMPLEX, 0.7, (255, 255, 255), 1)
        cv2.putText(img, str(round(section_acier_x, 2)), (200, 220), cv2.FONT_HERSHEY_COMPLEX, 0.7, (255, 255, 255), 1)
        cv2.putText(img, "Y-steel :", (360, 220), cv2.FONT_HERSHEY_COMPLEX, 0.7, (255, 255, 255), 1)
        cv2.putText(img, str(round(section_acier_y, 2)), (520, 220), cv2.FONT_HERSHEY_COMPLEX, 0.7, (255, 255, 255), 1)
    cv2.putText(img, "spacing_min_x :", (50, 250), cv2.FONT_HERSHEY_COMPLEX, 0.7, (255, 255, 255), 1)
    cv2.putText(img, str(round(espacement_x, 2)), (250, 250), cv2.FONT_HERSHEY_COMPLEX, 0.7, (255, 255, 255), 1)
    cv2.putText(img, "spacing_min_y :", (400, 250), cv2.FONT_HERSHEY_COMPLEX, 0.7, (255, 255, 255), 1)
    cv2.putText(img, str(round(espacement_y, 2)), (600, 250), cv2.FONT_HERSHEY_COMPLEX, 0.7, (255, 255, 255), 1)


    def affiche():
        st.image(image=img, output_format="auto")
        print(espacement_x)
        print(espacement_y)


    resultat = st.sidebar.button("result", on_click = affiche)
