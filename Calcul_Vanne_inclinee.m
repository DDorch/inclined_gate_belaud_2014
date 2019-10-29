%-------------------------------------------------------------------------------------------------------
% Programme de calcul d'une vanne inclinee
% Auteur: G. Belaud (UMR G-eau) - février 2016
% Mise à jour 01 Juin 2017: calcul pour vanne complètement ouverte, avec loi universelle de perte de charge
% Mise à jour 28 Aoùt 2018: ajout fonction, révision test free flow
% Mise à jour 31 Octobre 2018: calcul de Cc de façon tabulaire
%-------------------------------------------------------------------------------------------------------

test_vanne("data_test", 1);
test_vanne("test45deg", 1);
test_vanne("data_test", 2);
test_vanne("test45deg", 2);

% Notations:
% W=ouverture
% 0: section amont ; H: charge, h hauteur d'eau
% h1=hauteur d'eau au niveau de la section contractée
% s=h1/H0
% a=W/H0

% Fonctions externes à utiliser:
% CcFree : calcule l'angle pour une vanne dénoyée à faible ouverture, en
% fonction de l'angle

function test_vanne(sInput, Modele)
    %-------------------------------------------------------------------------------------------------------
    % Initialisation, constantes, parametres de calcul
    %-------------------------------------------------------------------------------------------------------

    %%%%%geometrie_vanne; % definition de la geometrie de la vanne, autres constantes
    Cc0 = 0.61;
    corr = 0; %coefficient correcteur (pertes de charge, etc.) - Ref: Belaud et al., (2009), JHE [coefficent k]
    g = 9.81; %constante gravitationnelle

    %-------------------------------------------------------------------------------------------------------
    % Niveaux d'eau, ouvertures et positions. Valeurs en mètres pour les
    % longueurs, en degres pour les angles
    %-------------------------------------------------------------------------------------------------------
    %X=load('data_test.txt','-ascii');
    X = load(strcat(sInput, '.txt'), '-ascii');
    %X=load('test45degFree.txt','-ascii');
    Zam = X(:, 1); %cote amont - reference=crete ouvrage
    Zav = X(:, 2); %cote aval - reference=crete ouvrage
    Wvanne = X(:, 3); %ouverture
    Avanne = X(:, 4); %angle en degres - angle entre l'horizontale (orientée vers l'amont) et la vanne (90=verticale, <90 ouverte vers l'aval)
    Nmes = length(Zam);

    % donnees de la vanne fixes
    pelle_amont = 0; % hauteur entre la crete (=cote 0) et le fond du canal à l'amont. Pourra etre mis comme vecteur si besoin
    pelle_aval = 0; % hauteur de la marche côté aval de la vanne (decrochement vers le bas)
    B = 1; % largeur de la vanne
    B2 = 1; % largeur du canal aval
    DX = 10; % distance en m entre les mesures amont et aval - utile uniquement pour vanne compl?tement ouverte pour calculer une perte de charge
    Cf = 0.02; % coefficient de frottement pour le calcul de cette perte de charge

    %-------------------------------------------------------------------------------------------------------
    N1 = 1; %indice de depart de la ligne de donnees a calculer
    N2 = Nmes; %indice de fin de la ligne de donnees a calculer

    %-------------------------------------------------------------------------------------------------------
    % Sauvegarde des résultats dans un fichier texte
    %-------------------------------------------------------------------------------------------------------
    sModele = ["formule", "tabul"];
    u1 = fopen(strcat('results_calcul_Q_v20181031_', sInput, '_', sModele(Modele), '.txt'), 'w');
    fprintf(u1, '%s\n', 'Results of gate calculation, according to the method published in Belaud et al, 2014, with tabulated Cc');
    fprintf(u1, '%s\n', '    h0     h2      W angle Dischar. R     Cd     Cc     h1  dQ/dW dQ/dh0 dQ/dh2 nbIter -- R: regime noye=3, partiellement noye=2, denoye=1 ');

    %-------------------------------------------------------------------------------------------------------
    % Boucle sur les differentes donnees
    %-------------------------------------------------------------------------------------------------------
    for k = N1:N2
        h0 = Zam(k) + pelle_amont; % hauteur amont
        z2 = Zav(k); %cote aval=hauteur par rapport à la crete
        W = Wvanne(k); % ouverture
        angle = Avanne(k); % angle de la vanne

        [Q, h1, etatVanne, Cd, Cc_app, dQdW, dQdhm, dQdhv, nbIter] = CalculQ(h0, z2, pelle_aval, W, angle, Cc0, B, B2, corr, Modele);

        fprintf(u1, '%6.3f %6.3f %6.3f %5.1f %8.4f %1i %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %3i\n', h0, z2 + pelle_aval, W, angle, Q, ...
            etatVanne, Cd, Cc_app, h1, dQdW, dQdhm, dQdhv, nbIter);
    end

    fclose(u1);
end

function [Q, h1, etatVanne, Cd, Cc_app, dQdW, dQdhm, dQdhv, nbIter] = CalculQ(h0, z2, pelle_aval, W, angle, Cc0, B, B2, corr, Modele)

    Imax = 20; % nombre d'iterations de boucles maximal
    epsQ1 = 0.001; % precision relative sur le debit pour le calcul itératif

    g = 9.81; %constante gravitationnelle

    h2 = z2 + pelle_aval; % cote aval par rapport au radier

    % ----------------Vanne inclinee-------------------
    Qref = Cc0 * W * B * sqrt(2 * g * h0); %Valeur de refernce pour le debit (pour tests d'arret)
    dQmax = max(0.001, Qref) * epsQ1; % precision attendue sur Q
    Q1 = Qref;

    if W < (h0 * 0.99)
        %------------------------------------------------------------
        % Vanne touchant ou rentrant dans l'eau
        %------------------------------------------------------------

        % Initialisations
        Free = 0;
        H0 = h0; % Point de depart initialisation: charge=hauteur
        Cc_F = Cc(W / H0, 0, angle, Modele); %appel de la fonction approchee du calcul de Cc

        h1_F = Cc_F * W; %hauteur definie par le jet issu de la vanne d?noyee

        if z2 <= h1_F
            % L'ecoulement est denoy?: pas d'influence aval
            Free = 1;
            h1 = h1_F;
        else
            Free = 0; % on n'est pas sur que l'ecoulement soit libre
            h1 = 0.5 * (z2 + h1_F); % valeur a priori, comme moyenne entre les 2 niveaux h1_F et z2
        end

        % Boucle itérative en supposant H0=h0
        i = 0;
        dQ1 = 1000;

        while (i < Imax) & (abs(dQ1) > dQmax)
            i = i + 1; % compteur de boucle
            % valeurs initiales estimees
            a = W / H0;
            %Test si le niveau aval est

            s = h1 / H0;
            % Hauteur de la veine contractee
            Cc1 = Cc(a, s, angle, Modele); %appel de la fonction approchee du calcul de Cc - a verifier pour angles>90
            h3 = Cc1 * W; % hauteur contractee avec a, s et H0 estimes

            % Debit estime par Bernoulli entre amont et section contractee, avec
            % correction donnee par eps
            Qit = sqrt(2 .* g * (h0 - h1) / (1 ./ h3^2 - 1 ./ h0^2)) * B / sqrt(1 + corr);
            h1it = h3;

            %Fin du calcul si l'écoulement est denoyé

            if Free == 0
                %On est dans le cas incertain noyé ou denoye, mais z2=h2-pelle_aval >h1_F
                % Tests de quantite de mouvement entre section contractee et aval: M2=Impulsion aval
                % Attention, en cas de pelle aval non nulle, M2 tient compte du
                % fait que la force sur la pelle aval s'annule
                M2 = g * B2 * z2^2 + 2 .* Qit^2 / (h2 * B2);
                % M3=impulsion dans la veine contractee, calculée avec la vitesse dans le
                % jet et la force de pression hydrostatique de la colonne de hauteur h1
                M3 = g * B2 * h1^2 + 2 * Qit^2 / (h3 * B); % expression qu'on pourra faire evoluer pour tenir compte de la hauteur d'eau
                % sur les cotes, dans le cas ou la vanne est moins large que l'ecoulement (B2>B)

                h3m = Cc(a, 0, angle, Modele) * W; % Hauteur contractee en ecoulement libre
                %impulsion aval inf?rieure à l'impulsion amont
                if (M2 < M3)
                    h1it = h3m; % On fait l'hypothese de free flow car le ressaut va alors se faire apres une courbe de remous F3
                    Free = 1;
                else
                    % Impulsion aval superieure à l'impulsion amont
                    % M3 insuffisant: il faut rajouter de la force de pression (h1 doit être plus élevé) pour equilibrer M2
                    % hauteur permettant de retrouver l'impulsion aval, en supposant Q
                    % inchangé:
                    dM = (M2 - 2 .* Qit^2 / (h3 * B)) / (g * B);

                    if (dM < h3m^2)% cas qui ne devrait pas arriver
                        h1it = h3m; % Niveau denoye - niveau physiquement minimum
                        Free = 1;
                    else
                        h1it = sqrt(dM); % noye
                    end

                end

            end

            % Actualisation
            %dh1=(h1it-h1)/h1;
            %h1=h1it;
            dQ1 = (Qit - Q1) / Qref; %test d'arret sur Q
            Q1 = Qit;

            Qt(i) = Qit; % valeurs temporaires pour tester la convergence
            h1t(i) = h1;
            H0 = h0 + Qit^2 / (2 .* g * B^2 * h0^2); %mise a jour de la charge amont

            %------------------------------
            % Fin des iterations
            %------------------------------
        end

        % Fin du calcul pour l'ouverture w
        %/ Sauvegarde des valeurs obtenues
        a = W / H0;
        s = h1 / H0;
        nbIter = i; %nombre d'iterations
        Cd = Qit / (B * W * sqrt(2 * g * h0)); % coefficient de debit
        Cc_app = Cc1; % coefficient de contraction sous la vanne
        h1c = h1; % hauteur d'eau au centre, juste apres la vanneQ=Qit;
        Q = Qit;

        %----------------------
        % Calcul des dérivées
        %----------------------
        % Etat de la vanne (noye=3, partiellement noye=2, denoye=1)
        F2 = (Qit / B)^2 / (g * h0^3);
        dQdhm = Qit / (2 * (H0 - h1)) * (1 - F2);
        dQdW = Qit / W;

        if s > a
            etatVanne = 3;
            dQdhv = -dQdhm;
        elseif s > Cc1
            etatVanne = 2;
            dQdhv = -dQdhm;
        else
            etatVanne = 1;
            dQdhv = 0;
        end

    else
        %-----------------------------------------------
        % vanne completement ouverte
        % Cette procedure etait implementee pour la VPT
        % Voir comment le gerer, avec le minimum de donnees: celles qui permettent
        % d'evaluer la perte de charge dans une vanne ouverte (loi de
        % retrecissement-elargissement par exemple, ou Manning si pas de
        % retrecissement.
        % Si on peut recuperer le debit au pas de temps anterieur, on peut le
        % recuperer, sinon mettre recuperer les parametres des sections du biefs,
        % et calculer loi de retrecissement-elergissement selon Lencastre
        % Ecrire algo qui determine Q tq DeltaH=perte de charge lineaire sur DX
        % metres +singuliere liee a elargissement+retrecissement
        % Dans ce code il n'y a pas les param?tres n?cessaires ni toutes les fonctions
        %-----------------------------------------------
        etatVanne = 4;
        %    if k>1
        %        %QVPT=QVPT(k-1);
        %        hv=(h0+h2)/2.;
        %        S=trapeze(Lfond,Bslope,hv);
        %        Q=Ch*S^1.5*sqrt((h0-h2)/((Lfond+2*sqrt(1+Bslope^2)*hv)*DX));
        %    else
        %        Q=0;
        %    end

        % Hypothèse: loi de frottement universelle, donnant la vitesse
        dH = Ham - Hav;
        hmoy = 0.5 * (Ham + Hav);
        Smoy = hmoy * B;
        R = Smoy / (B + 2 .* hmoy);
        U = sqrt(8 .* R * g * dH / (Cf * DX));
        Q = U * Smoy;
        dQdW = 0.;
        dQdhm = Smoy * sqrt(2 .* R * g / (Cf * DX * dH));
        dQdhv = -dQdhm;

    end
end
