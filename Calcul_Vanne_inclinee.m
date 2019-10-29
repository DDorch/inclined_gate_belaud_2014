%-------------------------------------------------------------------------------------------------------
% Programme de calcul d'une vanne inclinee
% Auteur: G. Belaud (UMR G-eau) - février 2016
% Mise à jour 01 Juin 2017: calcul pour vanne complètement ouverte, avec loi universelle de perte de charge
%-------------------------------------------------------------------------------------------------------

% Notations:
% W=ouverture
% 0: section amont ; H: charge, h hauteur d'eau
% h1=hauteur d'eau au niveau de la section contractée
% s=h1/H0
% a=W/H0

%-------------------------------------------------------------------------------------------------------
% Initialisation, constantes, parametres de calcul
%-------------------------------------------------------------------------------------------------------

%%%%%geometrie_vanne; % definition de la geometrie de la vanne, autres constantes
Imax = 50; % nombre d'iterations de boucles maximal
epsh1 = 0.001; % precision relative sur le niveau aval pour le calcul itÃ©ratif
Cc0 = 0.61;
corr = 0; %coefficient correcteur (pertes de charge, etc.)
g = 9.81; %constante gravitationnelle

%-------------------------------------------------------------------------------------------------------
% Niveaux d'eau, ouvertures et positions. Valeurs en mï¿½tres pour les
% longueurs, en degres pour les angles
%-------------------------------------------------------------------------------------------------------
X = load('data_test.txt', '-ascii');
Ham = X(:, 1); %hauteur amont
Hav = X(:, 2); %hauteur aval
Wvanne = X(:, 3); %ouverture
Avanne = X(:, 4); %angle en degres - angle entre l'horizontale et la vanne (90=verticale, <90 ouverte vers l'aval)
Nmes = length(Ham);

% donnees de la vanne fixes
cote_radier = 0; % pourra etre mis comme vecteur si besoin
B = 1;
B2 = 1;
DX = 10; % distance entre les mesures amont et aval - utile uniquement pour vanne complètement ouverte
Cf = 0.02; % coefficient de frottement

%-------------------------------------------------------------------------------------------------------
N1 = 1;
N2 = Nmes;

for k = N1:N2
    %-------------------------------------------------------------------------------------------------------
    % Boucle sur les differentes mesures
    %-------------------------------------------------------------------------------------------------------
    h0 = Ham(k) - cote_radier; % hauteur amont
    h2 = Hav(k) - cote_radier; % hauteur aval

    % ----------------Vanne inclinee-------------------

    W = Wvanne(k); % ouverture
    angle = Avanne(k); % angle de la vanne

    if W < (h0 * 0.99)
        %------------------------------------------------------------
        % Vanne touchant ou rentrant dans l'eau
        %------------------------------------------------------------

        % Initialisations
        H0 = h0; % Point de depart initialisation: charge=hauteur
        h1 = 0.5 * (h2 + Cc0 * W); % valeur a priori, comme moyenne entre un niveau denoye et completement noye
        h1 = h2;

        % Boucle itérative
        i = 0;
        dh1 = 1000;

        while (i < Imax) & (abs(dh1) > epsh1)
            i = i + 1; % compteur de boucle
            % valeurs initiales estimees
            a = W / H0;
            s = h1 / H0;
            % Hauteur de la veine contractee
            Cc1 = Cc(a, s, angle); %appel de la fonction approchee du calcul de Cc - a verifier pour angles>90
            h3 = Cc1 * W; % hauteur contractee avec a, s et H0 estimes
            h3m = Cc(a, 0, angle) * W; % Hauteur contractee en ecoulement libre

            % Debit estime par Bernoulli entre amont et section contractee, avec
            % correction donnee par eps
            Qit = sqrt(2 .* g * (h0 - h1) / (1 ./ h3^2 - 1 ./ h0^2)) * B / sqrt(1 .+ corr);

            % Tests de quantite de mouvement entre section contractee et aval:
            %Impulsion aval
            M2 = g * B2 * h2^2 + 2 .* Qit^2 / (h2 * B2);
            % impulsion dans la veine contractee
            M3 = g * B2 * h1^2 + 2 * Qit^2 / (h3 * B); % expression qu'on pourra faire evoluer pour tenir compte de la hauteur d'eau
            % sur les cotes, dans le cas ou la vanne est moins large que
            % l'ecoulement (B2>B)

            %impulsion aval inférieure à l'impulsion amont
            if (M2 < M3)
                h1it = h3m; % On fait l'hypothese de free flow car le ressaut va alors se faire apres une courbe de remous F3
            else
                % M3 insuffisant: il faut rajouter de la force de pression pour equilibrer M2
                % hauteur permettant de retrouver l'impulsion aval
                dM = (M2 - 2 .* Qit^2 / (h3 * B)) / (g * B);

                if (dM < h3m^2)
                    h1it = h3m; % Niveau denoye - niveau physiquement minimum
                else
                    h1it = sqrt(dM); % noye
                end

            end

            % Actualisation
            dh1 = (h1it - h1) / h1;
            h1 = h1it;
            Qt(i) = Qit; % valeurs temporaires pour tester la convergence
            h1t(i) = h1;
            H0 = h0 + Qit^2 / (2 .* g * B^2 * h0^2); %mise a jour de la charge amont

            % Fin du calcul pour l'ouverture w(k)
        end

        %/ Sauvegarde des valeurs obtenues
        a = W / H0;
        s = h1 / H0;
        NbIter(k) = i; %nombre d'iterations
        Cd(k) = Qit / (B * W * sqrt(2 * g * h0)); % coefficient de debit
        Cc_app(k) = Cc1; % coefficient de contraction sous la vanne
        h1c(k) = h1; % hauteur d'eau au centre, juste apres la vanneQ(k)=Qit;
        Q(k) = Qit;

        %----------------------
        % Calcul des dérivées
        %----------------------
        % Etat de la vanne (noye=3, partiellement noye=2, denoye=1)
        F2 = (Qit / B)^2 / (g * h0^3);
        dQdhm = Qit / (2 * (H0 - h1)) * (1 - F2);
        dQdW = Qit / W;

        if s > a
            etatVanne(k) = 3;
            dQdhv = -dQdhm;
        elseif s > Cc1
            etatVanne(k) = 2;
            dQdhv = -dQdhm;
        else
            etatVanne(k) = 1;
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
        % Dans ce code il n'y a pas les paramètres nécessaires ni toutes les fonctions
        %-----------------------------------------------
        etatVanne(k) = 4;
        %    if k>1
        %        %QVPT(k)=QVPT(k-1);
        %        hv=(h0+h2)/2.;
        %        S=trapeze(Lfond,Bslope,hv);
        %        Q(k)=Ch*S^1.5*sqrt((h0-h2)/((Lfond+2*sqrt(1+Bslope^2)*hv)*DX));
        %    else
        %        Q(k)=0;
        %    end

        % Hypothèse: loi de frottement universelle, donnant la vitesse
        dH = Ham - Hav;
        hmoy = 0.5 * (Ham + Hav);
        Smoy = hmoy * B;
        R = Smoy / (B + 2 .* hmoy);
        U = sqrt(8 .* R * g * dH / (Cf * DX));
        Q(k) = U * Smoy;
        dQdW = 0.;
        dQdhm = Smoy * sqrt(2 .* R * g / (Cf * DX * dH));
        dQdhv = -dQdhm;
    end

end
