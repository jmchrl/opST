#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2022 Jean-Michel CHEREL
#
# This file is part of pyBar.
#    pyBar is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 3 of the License, or
#    (at your option) any later version.
#
#    pyBar is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with pyBar; if not, write to the Free Software
#    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

import classRdm as pB
import function
import math
import numpy
import xml.etree.ElementTree as ET

CSECTIONS = "c:/users/jmc/pyBar/library/section.xml"
CMTX = "c:/users/jmc/pyBar/library/material.xml"

def efforts_barre(chargement, barre, l, rdm):
    """Retourne les efforts le long d'une barre (N, V, Mf) sous forme de tableau pour 1 chargement"""

    charFp = chargement.charBarFp.get(barre, {})
    charQu = chargement.charBarQu.get(barre, {})
    charTri = chargement.charBarTri.get(barre, {})
    tenseur = numpy.zeros((3, 11), dtype = float)
    x = 0.00
    colonne = 0
    while x <= 1:
        tenseur[0, colonne] = rdm.NormalPoint(chargement, barre, x*l, charQu, charFp, charTri, rdm.conv)
        tenseur[1, colonne] = rdm.TranchantPoint(chargement, barre, x*l, charQu, charFp, charTri, rdm.conv)
        tenseur[2, colonne] = rdm.MomentPoint(chargement, barre, x*l, charQu, charFp, charTri, rdm.conv)
        x = x + 0.10
        colonne = colonne + 1
    return tenseur

def verif_barre(barre, l, section_barre, mtx_barre, efforts_combines, lgfz, lgfy, courbef, ldfz, v_deversement):
    """Verification des barres suivants les eurocodes"""
    
    barre = BarreEC3(barre, l, section_barre, mtx_barre, efforts_combines, lgfz, lgfy, courbef, ldfz, v_deversement)
    return barre

class BarreEC3():
    """Verification d une barre suivant l eurocode 3"""

    def __init__(self, barre, l, section_barre, mtx_barre, efforts_combines, lgfz, lgfy, courbef, ldfz, v_deversement):
        
        self.nom = barre
        self.longueur = l
        self.nom_section = section_barre
        self.nom_mtx = mtx_barre
        self.efforts = efforts_combines
        self.young = None
        self.fy = None
        self.__recherche_caracteristiques_mtx()
        self.s = None
        self.igz = None
        self.h = None
        self.b = None
        self.vz = None
        self.igy = None
        self.vy = None
        self.ifz = None
        self.__recherche_caracteristiques_section()
        self.rayon_giration_z = math.sqrt(self.igz/self.s)
        self.rayon_giration_y = math.sqrt(self.igy/self.s)
        self.lgfz = lgfz
        self.lgfy = lgfy
        self.courbef = courbef
        self.ldfz = ldfz
        self.alpha = None
        self.__facteur_imperfection()
        self.v_deversement = v_deversement
        self.gammaM0 = 1
        self.gammaM1 = 1
        
        # Lancement de la routine de verification de la barre à l ELU
        
        self.resultatELU = self.__routine_verif_barre_ELU()
        print(self.resultatELU)
    
    def __recherche_caracteristiques_mtx(self):
        
        xml_mtx = ET.parse(CMTX)
        strxml_mtx = xml_mtx.getroot()
        for enfant in strxml_mtx:
            if enfant.tag == "group":
                for petit_enfant in enfant:
                    if petit_enfant.get('id') == self.nom_mtx:
                        self.young = float(petit_enfant.get('young'))
                        self.fy = float(petit_enfant.get('fy'))
            else:
                if enfant.get('id') == self.nom_mtx:
                    self.young = float(enfant.get('young'))
                    self.fy = float(enfant.get('fy'))
    
    def __recherche_caracteristiques_section(self):
        
        xml_section = ET.parse(CSECTIONS)
        strxml_section = xml_section.getroot()
        for enfant in strxml_section:
            if enfant.tag == "group":
                for petit_enfant in enfant:
                    if petit_enfant.get('id') == self.nom_section:
                        self.s = float(petit_enfant.get('s'))
                        self.h = float(petit_enfant.get('h'))
                        self.b = float(petit_enfant.get('b'))
                        # Attention a la convention d axe dans pyBar, igy dans pyBar et section.xml est igz pour la suite des verifications
                        self.igz = float(petit_enfant.get('igy'))
                        # Attention a la convention d axe dans pyBar, vz dans pyBar et section.xml est vy pour la suite des verifications
                        self.vy = float(petit_enfant.get('vz'))
                        # Attention a la convention d axe dans pyBar, igz dans pyBar et section.xml est igy pour la suite des verifications
                        self.igy = float(petit_enfant.get('igz'))
                        # Attention a la convention d axe dans pyBar, v dans pyBar et section.xml est vz pour la suite des verifications
                        self.vz = float(petit_enfant.get('v'))
                        # Attention a la convention d axe dans pyBar, ify dans pyBar et section.xml est ifz pour la suite des verifications
                        self.ifz = float(petit_enfant.get('ify'))
            else:
                if enfant.get('id') == self.nom_section:
                    self.s = float(petit_enfant.get('s'))
                    self.h = float(petit_enfant.get('h'))
                    self.b = float(petit_enfant.get('b'))
                    # Attention a la convention d axe dans pyBar, igy dans pyBar et section.xml est igz pour la suite des verifications
                    self.igz = float(petit_enfant.get('igy'))
                    # Attention a la convention d axe dans pyBar, vz dans pyBar et section.xml est vy pour la suite des verifications
                    self.vy = float(petit_enfant.get('vz'))
                    # Attention a la convention d axe dans pyBar, igz dans pyBar et section.xml est igy pour la suite des verifications
                    self.igy = float(petit_enfant.get('igz'))
                    # Attention a la convention d axe dans pyBar, v dans pyBar et section.xml est vz pour la suite des verifications
                    self.vz = float(petit_enfant.get('v'))
                    # Attention a la convention d axe dans pyBar, ify dans pyBar et section.xml est ifz pour la suite des verifications
                    self.ifz = float(petit_enfant.get('ify'))

    def __facteur_imperfection(self):
        
        if self.courbef == 'a0':
            self.alpha = 0.13
        if self.courbef == 'a':
            self.alpha = 0.21
        if self.courbef == 'b':
            self.alpha = 0.34
        if self.courbef == 'c':
            self.alpha = 0.49
        if self.courbef == 'd':
            self.alpha = 0.76

    def __routine_verif_barre_ELU(self):
        
        resultat = {}
        resultat['ratio_section'] = 0
        resultat['ratio_flambement'] = 0
        for cle in list(self.efforts.keys()):
            if cle[:3] == "ELU":
                point = 0
                indice = 0
                while point <= 1:
                    if self.efforts[cle][0,indice] != 0.0:
                        if self.efforts[cle][2,indice] == 0.0:
                            # Compression simple ou traction simple
                            ratio = self.verification_Traction_Compression(self.efforts[cle][0,indice])
                            if ratio > resultat['ratio_section']:
                                resultat['ratio_section'] = round(ratio,2)
                                resultat['point'] = indice
                                resultat['position'] = str(point) + '*l = ' + str(round(point*self.longueur,3))
                                resultat['combinaison'] = cle
                                resultat['Ned'] = round(self.efforts[cle][0,indice],2)
                                resultat['Ved'] = round(self.efforts[cle][1,indice],2)
                                resultat['Med'] = round(self.efforts[cle][2,indice],2)
                            if self.efforts[cle][0,indice] < 0.0:
                                # Flambement
                                ratio = self.flambement(self.efforts[cle][0,indice])
                                if ratio > resultat['ratio_flambement']:
                                    resultat['ratio_flambement'] = round(ratio,2)
                                    resultat['point'] = indice
                                    resultat['position'] = str(point) + '*l = ' + str(round(point*self.longueur,3))
                                    resultat['combinaison'] = cle
                                    resultat['Ned'] = round(self.efforts[cle][0,indice],2)
                                    resultat['Ved'] = round(self.efforts[cle][1,indice],2)
                                    resultat['Med'] = round(self.efforts[cle][2,indice],2)
                        else:
                            # Flexion composee
                            ratio = self.verification_flexion_composee(self.efforts[cle][2,indice],self.efforts[cle][0,indice])
                            if ratio > resultat['ratio_section']:
                                resultat['ratio_section'] = round(ratio,2)
                                resultat['deversement'] = self.deversement(self.efforts[cle][2,indice], self.efforts[cle])
                                resultat['point'] = indice
                                resultat['position'] = str(point) + '*l = ' + str(round(point*self.longueur,3))
                                resultat['combinaison'] = cle
                                resultat['Ned'] = round(self.efforts[cle][0,indice],2)
                                resultat['Ved'] = round(self.efforts[cle][1,indice],2)
                                resultat['Med'] = round(self.efforts[cle][2,indice],2)
                    if self.efforts[cle][2,indice] != 0.0:
                        if self.efforts[cle][0,indice] == 0.0:
                            # Flexion simple
                            ratio = self.verification_flexion_simple(self.efforts[cle][2,indice])
                            if ratio > resultat['ratio_section']:
                                resultat['ratio_section'] = round(ratio,2)
                                resultat['deversement'] = self.deversement(self.efforts[cle][2,indice], self.efforts[cle])
                                resultat['point'] = indice
                                resultat['position'] = str(point) + '*l = ' + str(round(point*self.longueur,3))
                                resultat['combinaison'] = cle
                                resultat['Ned'] = round(self.efforts[cle][0,indice],2)
                                resultat['Ved'] = round(self.efforts[cle][1,indice],2)
                                resultat['Med'] = round(self.efforts[cle][2,indice],2)
                    point = point + 0.1
                    indice = indice + 1
        return resultat

    def verification_Traction_Compression(self, Ned):
        """verification en traction ou en compression, section de classe 3"""
        
        NplRd = (self.s*self.fy)/self.gammaM0
        ratio = abs(Ned)/NplRd
        return ratio
    
    def verification_flexion_simple(self, Med):
        """verification en flexion simple, section de classe 3"""
        
        MelRd = ((self.igy/self.vz)*self.fy)/self.gammaM0
        ratio = abs(Med)/MelRd
        return ratio

    def verification_flexion_composee(self, Med, Ned):
        """verification en flexion composee, section de classe 3"""
        MelRd = ((self.igy/self.vz)*self.fy)/self.gammaM0
        ratio = abs(Ned)/(self.s*self.fy) + abs(Med)/MelRd
        return ratio

    def elancement_reduit(self):
        """Calcul de l elancement reduit de la barre"""
        
        lambda_1 = math.pi*math.sqrt(self.young/self.fy)
        lambda_barre_y = (self.lgfy/self.rayon_giration_y)*(1/lambda_1)
        lambda_barre_z = (self.lgfz/self.rayon_giration_z)*(1/lambda_1)
        return max(lambda_barre_y, lambda_barre_z)

    def flambement(self, Ned):
        """Verification de la resistance au flambement d une barre uniformement comprimee, section de classe 3"""
        
        lambda_barre = self.elancement_reduit()
        if lambda_barre <= 0.2:
            return 0
        else:
            phi = 0.5*(1+self.alpha*(lambda_barre-0.2)+lambda_barre**2)
            facteur_reduction_chi = min(1, 1/(phi+math.sqrt(phi**2-lambda_barre**2)))
            Nbrd = (facteur_reduction_chi*self.s*self.fy)/self.gammaM1
            ratio = abs(Ned)/self.gammaM1
            return ratio

    def elancement_limite_semelle_comprimee(self):
        """Calcul de l elencement limite de semelle comprimee equivalente, cf. 6.3.2.4 EC3"""
        
        if self.nom_section[0] == "I" or "H":
            elancement = 0.2+0.1*(self.b/self.h)
        if self.nom_section[:3] == "PRS":
            elancement = 0.3*(self.b/self.h)
        else:
            elancement = 0.2
        return elancement

    def deversement(self, Med, efforts):
        """Verification de la sensibilite au deversement suivant methode simplifiee, cf. 6.3.2.4 EC3"""
        
        lambda_barre_c0 = self.elancement_limite_semelle_comprimee()
        
        # Calcul de Mcrd
        
        if self.v_deversement == "v_max":
            Wy = self.igy/self.vy
        if self.v_deversement == "v_min":
            Wy = self.igy/(self.h-self.vy)
        else:
            Wy = self.igy/(self.h/2)
        
        Mcrd = (Wy*self.fy)/self.gammaM1

        # Calcul de kc
        
        if efforts[2,0] >0.001 and efforts[2,10] >0.001 and efforts[2,5] >0.001:
            if efforts[2,0] == efforts[2,10] and efforts[2,0] == efforts[2,5]:
                kc = 1
            else:
                psi = 1/(efforts[2,0]/efforts[2,10])
                kc = 1/(1.33-0.33*psi)
        if efforts[2,0] <-0.001 and efforts[2,10] <-0.001 and efforts[2,5] <-0.001:
            if efforts[2,0] == efforts[2,10] and efforts[2,0] == efforts[2,5]:
                kc = 1
            else:
                psi = 1/(efforts[2,0]/efforts[2,10])
                kc = 1/(1.33-0.33*psi)
        if efforts[2,0] >0.001 and efforts[2,10] <-0.001:
            if efforts[2,0] - efforts[2,5] == efforts[2,5] - efforts[2,10]:
                psi = 1/(efforts[2,0]/efforts[2,10])
                kc = 1/(1.33-0.33*psi)
        
        if efforts[2,0] <0.001 and efforts[2,0] > -0.001 and efforts[2,10] <0.001 and efforts[2,10]:
            if efforts[2,0] - efforts[2,1] == efforts[2,1] - efforts[2,2]:
                kc = 0.86
            else:
                kc = 0.94
        
        if efforts[2,0] > 0.001 and efforts[2,10] >0.001 and efforts[2,5] < 0.001:
            if efforts[2,0] - efforts[2,1] == efforts[2,1] - efforts[2,2]:
                kc = 0.7
            else:
                kc = 0.90
        
        if efforts[2,0] <0.001 and efforts[2,0] > -0.001 and efforts[2,10] > 0.001 and efforts[2,5] < 0.001:
            if efforts[2,0] - efforts[2,1] == efforts[2,1] - efforts[2,2]:
                kc = 0.82
            else:
                kc = 0.91
        
        if efforts[2,0] <0.001 and efforts[2,0] > -0.001 and efforts[2,10] < 0.001 and efforts[2,5] > 0.001:
            if efforts[2,0] - efforts[2,1] == efforts[2,1] - efforts[2,2]:
                kc = 0.82
            else:
                kc = 0.91
        
        if efforts[2,10] <0.001 and efforts[2,10] > -0.001 and efforts[0,2] > 0.001 and efforts[2,10] < 0.001:
            if efforts[2,0] - efforts[2,1] == efforts[2,1] - efforts[2,2]:
                kc = 0.82
            else:
                kc = 0.91
        
        if efforts[2,10] <0.001 and efforts[2,10] > -0.001 and efforts[2,0] < 0.001 and efforts[2,5] > 0.001:
            if efforts[2,0] - efforts[2,1] == efforts[2,1] - efforts[2,2]:
                kc = 0.82
            else:
                kc = 0.91
        
        else:
            input('La valeur de kc n a pas ete trouvee, indiquer la valeu de kc\n cf. page 151 precis genie-civil kc :')
        
        # Calcul de lambda_1
        
        lambda_1 = math.pi*math.sqrt(self.young/self.fy)
        
        # Calcul de lambda barre f
        
        lambda_barre_f = (kc*self.ldfz)/(self.ifz*lambda_1)
        
        # Verification
        
        if lambda_barre_f <= lambda_barre_c0*(Mcrd/Med):
            return 'Pas de risque de deversement'
        else:
            return 'Attention risqe de deversement'
        
        

if __name__ == "__main__":
    file = input("chemin du fichier de donnee pyBar : ")
    xml = ET.parse(file)
    strxml = xml.getroot()
    xml_types_barres = ET.parse('opSTtypesBarres.dat')
    strxml_types_barres = xml_types_barres.getroot()
                
    structure = pB.Structure(xml, file)
    rdm = pB.R_Structure(structure)

    efforts_combines = {}
    tenseur_combi = numpy.zeros((3, 11), dtype = float)

    for barre in rdm.struct.GetBars():
        for nom_combi, coefficients in rdm.GetCombi().items():
            for nom_cas, coef_cas in coefficients.items():
                for numero_cas in range(len(rdm.GetCasCharge())):
                    if rdm.GetCharNameByNumber(numero_cas) == nom_cas:
                        tenseur_cas = efforts_barre(rdm.GetCharByNumber(numero_cas), barre, rdm.struct.Lengths[barre], rdm)
                tenseur_combi = tenseur_combi + coef_cas * tenseur_cas
            efforts_combines[nom_combi] = tenseur_combi
            tenseur_combi = numpy.zeros((3, 11), dtype = float)
        
        # recherche du nom de la section et du nom du materiaux de la barre
        for elem in strxml.findall('elem'):
            if elem.get('id') == "geo":
                for section in elem.findall('barre'):
                    if section.get('id') == "*":
                        section_barre = section.get('profil')
                    else:
                        barres_section = section.get('id').split(',')
                        for barre_section in barres_section:
                            if barre_section == barre:
                                section_barre = section.get('profil')
            if elem.get('id') == "material":
                for material in elem.findall('barre'):
                    if material.get('id') == "*":
                        mtx_barre = material.get('profil')
                    else:
                        barres_mtx = material.get('id').split(',')
                        for barre_mtx in barres_mtx:
                            if barre_mtx == barre:
                                mtx_barre = material.get('profil')
        
        # recherche des longueurs de flambement et type de courbe de flambement
        for type_de_barre in strxml_types_barres:
            if type_de_barre.get('barres') == "*":
                if type_de_barre.get('lgfz')[:4] == 'coef':
                    lgfz = float(type_de_barre.get('lgfz')[4:])* rdm.struct.Lengths[barre]
                else:
                    lgfz = float(type_de_barre.get('lgfz')[4:])
                if type_de_barre.get('lgfy')[:4] == 'coef':
                    lgfy = float(type_de_barre.get('lgfy')[4:])* rdm.struct.Lengths[barre]
                else:
                    lgfy = float(type_de_barre.get('lgfy')[4:])
                courbef = type_de_barre.get('courbef')
                if type_de_barre.get('ldfz')[:4] == 'coef':
                    ldfz = float(type_de_barre.get('ldfz')[4:])* rdm.struct.Lengths[barre]
                else:
                    ldfz = float(type_de_barre.get('ldfz')[4:])
                v_deversement = type_de_barre.get('v_deversement')
            else:
                type_barres = type_de_barre.get('barres').split(',')
                for type_barre in type_barres:
                    if type_barre == barre:
                        if type_de_barre.get('lgfz')[:4] == 'coef':
                            lgfz = float(type_de_barre.get('lgfz')[4:])* rdm.struct.Lengths[barre]
                        else:
                            lgfz = float(type_de_barre.get('lgfz')[4:])
                        if type_de_barre.get('lgfy')[:4] == 'coef':
                            lgfy = float(type_de_barre.get('lgfy')[4:])* rdm.struct.Lengths[barre]
                        else:
                            lgfy = float(type_de_barre.get('lgfy')[4:])
                        courbef = type_de_barre.get('courbef')
                        if type_de_barre.get('ldfz')[:4] == 'coef':
                            ldfz = float(type_de_barre.get('ldfz')[4:])* rdm.struct.Lengths[barre]
                        else:
                            ldfz = float(type_de_barre.get('ldfz')[4:])
                        v_deversement = type_de_barre.get('v_deversement')
        
        # Verification de la barre
        barre = verif_barre(barre, rdm.struct.Lengths[barre], section_barre, mtx_barre, efforts_combines, lgfz, lgfy, courbef, ldfz, v_deversement)
        
        # Edition de la note de calcul par barre
        template = open('latex/barre_EC3.tex', 'r')
        ndc_barre = open('ndc_barre.tex', 'w')
        note_barre = template.read()
        note_barre = note_barre.replace("NOMBARRE", str(barre.nom))
        note_barre = note_barre.replace("NUMPOINT", str(barre.resultatELU['point']))
        note_barre = note_barre.replace("POSITION", str(barre.resultatELU['position']))
        note_barre = note_barre.replace("COMBINAISON", str(barre.resultatELU['combinaison']))
        note_barre = note_barre.replace("MTX", str(barre.nom_mtx))
        note_barre = note_barre.replace("NOMSECTION", str(barre.nom_section))
        note_barre = note_barre.replace("S =", "S = " + str(barre.s))
        note_barre = note_barre.replace("I{g,z} =", "I{g,z} = " + str(barre.igz))
        note_barre = note_barre.replace("H =", "H = " + str(barre.h))
        note_barre = note_barre.replace("B =", "B = " + str(barre.b))
        note_barre = note_barre.replace("V{z} =", "V{z} = " + str(barre.vz))
        note_barre = note_barre.replace("I{g,y} =", "I{g,y} = " + str(barre.igy))
        note_barre = note_barre.replace("V{y} =", "V{y} = " + str(barre.igy))
        note_barre = note_barre.replace("I{f,z} =", "I{f,z} = " + str(barre.ifz))
        note_barre = note_barre.replace('EFFORTS', 'N{ed} = ' + str(barre.resultatELU['Ned']) + ' N\n\n' + \
                                                   'V{ed} = ' + str(barre.resultatELU['Ved']) + ' N\n\n' + \
                                                   'M{ed} = ' + str(barre.resultatELU['Med']) + ' N.m\n')
        ndc_barre.write(note_barre)