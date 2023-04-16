import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, roc_auc_score, roc_curve
import matplotlib.pyplot as plt
import seaborn as sns
    
def plot_roc(main_data, predictors:list, plot_title:str):
    # Load the dataset
    main_data.sex = main_data.sex.map({'M': 0, 'F': 1})

    # Define the features and target variable
    X = main_data[predictors]
    y = main_data['target_phenotype']

    # Split the dataset into training and testing sets
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.33, random_state=42, stratify = y)

    # Create a logistic regression model
    lr_model = LogisticRegression()

    # Train the model using the training data
    lr_model.fit(X_train, y_train)

    # Make predictions on the testing data
    y_pred = lr_model.predict(X_test)

    y_pred_proba = lr_model.predict_proba(X_test)[::,1]
    fpr, tpr, _ = roc_curve(y_test,  y_pred_proba)
    roc_auc = roc_auc_score(y_test, y_pred_proba)
    print("ROC AUC:", roc_auc)

    plt.plot(fpr,tpr, label='ROC curve (area = %0.5f)' % roc_auc)
    plt.title(plot_title)
    plt.plot([0, 1], [0, 1], 'k--') # Plot a diagonal line for reference
    plt.ylabel('True Positive Rate')
    plt.xlabel('False Positive Rate')
    plt.legend(loc="lower right")
    plt.show()
    
def plot_density(main_data, prs_column:str):
    # Create a figure and axis object
    fig, ax = plt.subplots()

    # Plot a histogram for each category
    sns.kdeplot(
        data = main_data, 
        x = prs_column, 
        hue = "target_phenotype", 
        fill = True, 
        common_norm = False,
        ax = ax)


    # Add axis labels and title
    ax.set_xlabel('PRS')
    ax.set_ylabel('Frequency')
    ax.set_title('PRS distribuition case vs control')