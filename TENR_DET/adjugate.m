%% Code for finding Adugate Matrix ..... using modified LU Decomposition .......

function  [Adjjacob]=adjugate(Jacob)

[L_Mat, U_Mat, P_Mat, Q_Mat]=lu(Jacob);              % LU Decomposition.....

[size_a,~]=size(Jacob);                       % Sizes....

diag_dvec= diag(L_Mat).*diag(U_Mat);          % Creatig D Matrix.....

Diag_Mat=sparse(diag(diag_dvec));             

New_LMat= zeros(size_a,size_a);

%% Creating new L Matrix

for i=1:size_a;
    New_LMat(:,i)=L_Mat(:,i)./L_Mat(i,i);
end

     New_LMat(eye(size(New_LMat))~=0)=0;
     
     New_LMat=New_LMat+eye(size_a,size_a);


%% Creating new U Matrix

New_UMat= zeros(size_a,size_a);

for i=1:size_a;
    New_UMat(i,:)=U_Mat(i,:)./U_Mat(i,i);
end
     New_UMat(eye(size(New_UMat))~=0)=0;
     
     New_UMat=New_UMat+eye(size_a,size_a);
%% Finding Adjugate Matrix for D ........


post_Diag_Mat=prod(diag_dvec)*eye(size_a,size_a);

Adjugate_Diagonal_L= diag(post_Diag_Mat)./diag_dvec;

Adjugate_DIAG=(det(Jacob)/det(Diag_Mat))*diag(Adjugate_Diagonal_L);


%% Adjugate Matrix Jacob......

UMMat=inv(New_UMat);
LLMat=inv(New_LMat);

Adjjacob=(UMMat*((Adjugate_DIAG*LLMat)))*inv(P_Mat);

%InvMatrix=(UMMat*((Diag_Mat*LLMat)))*P_Mat;
 
end
                                   