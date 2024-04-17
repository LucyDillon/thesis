from flask import Flask, render_template, flash, send_file, current_app, request
from flask_wtf import FlaskForm
from wtforms import FileField, SubmitField, SelectField, StringField, validators
from wtforms.validators import InputRequired
from snakemake import snakemake
from werkzeug.utils import secure_filename
import os
import uuid
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from email import encoders
import smtplib
import email_validator
from email.mime.base import MIMEBase
#from celery import Celery
import threading

# NEED TO ADD CORRECT URL: celery = Celery(__name__, broker='pyamqp://guest:guest@localhost//')

ALLOWED_EXTENSIONS = set(['fna', 'fasta'])

app = Flask(__name__)
app.config['SECRET_KEY'] = 'secretkey'
app.config['UPLOAD_FOLDER'] = 'static/files'
app.config['DOWNLOAD_FOLDER'] = 'static/files/results/'
# Set your email configurations
app.config['MAIL_SERVER'] = 'smtp.gmail.com'
app.config['MAIL_PORT'] = 587  # port for gmail
app.config['MAIL_USE_TLS'] = True
app.config['MAIL_USE_SSL'] = False
app.config['MAIL_USERNAME'] = 'amrphenopredict@gmail.com'
app.config['MAIL_PASSWORD'] = 'xxxxxxxxxxx' # Password removed from public view
app.config['MAIL_DEFAULT_SENDER'] = 'amrphenopredict@gmail.com'

class UploadFileForm(FlaskForm):
    file = FileField("File", validators=[InputRequired()])
    email = StringField("Email", validators=[validators.Email()])
    submit = SubmitField("Upload File")

def generate_unique_code():
    return str(uuid.uuid4().hex)

uploaded_files_by_ip = {}

def get_uploads_count_for_ip(ip):
    return uploaded_files_by_ip.get(ip, 0)

def increment_uploads_count_for_ip(ip):
    if ip in uploaded_files_by_ip:
        uploaded_files_by_ip[ip] += 1
    else:
        uploaded_files_by_ip[ip] = 1

def send_email(subject, recipients, body, attachments=None):
    msg = MIMEText(body)
    msg['Subject'] = subject
    msg['From'] = app.config['MAIL_DEFAULT_SENDER']
    msg['To'] = ', '.join(recipients)

    with smtplib.SMTP(app.config['MAIL_SERVER'], app.config['MAIL_PORT']) as server:
        server.starttls()
        server.login(app.config['MAIL_USERNAME'], app.config['MAIL_PASSWORD'])
        server.sendmail(app.config['MAIL_DEFAULT_SENDER'], recipients, msg.as_string())

snakemake_lock = threading.Lock()

@app.route('/', methods=['GET', 'POST'])
@app.route('/home', methods=['GET', 'POST'])
def home():
    form = UploadFileForm()
    result = None

    if form.validate_on_submit():
        ip = request.remote_addr

        # Check the upload count for this IP address
        uploads_count = uploaded_files_by_ip.get(ip, 0)
        max_uploads_per_ip = 20  # Set the maximum upload limit per IP address
        max_file_size = 13 * 1024 * 1024 # Max file size is set to 13 megabytes

        if uploads_count >= max_uploads_per_ip:
            return "Upload limit reached for this IP address", 403

        file = form.file.data
        if file:
            if len(file.read()) > max_file_size:
                return "File size exceeds the limit (13 MB)", 400  # Return an error response
            file.stream.seek(0)  # Reset the file stream position
            original_filename = secure_filename(file.filename)
            unique_code = generate_unique_code()
            filename_without_extension, file_extension = os.path.splitext(original_filename)
            unique_filename = f"{filename_without_extension}_{unique_code}.fna"
            file_path = os.path.join(app.config['UPLOAD_FOLDER'], unique_filename)
            file.save(file_path)

            selected_model = request.form.get('model')
            user_email = form.email.data
            if selected_model and user_email:
                with snakemake_lock:
                    if selected_model == 'DecisionTree':
                        snakefile = 'decision_tree.smk'
                    elif selected_model == 'CNN':
                        snakefile = 'cnn.smk'
                    elif selected_model == 'BothModels':
                        snakefile = 'both_models.smk'

                    snakemake(
                        snakefile=snakefile,
                        workdir="./",
                        config=dict(genome=filename_without_extension),
                        use_singularity=True,
                        cores=8
                    )

                result_filename = f"{filename_without_extension}_{unique_code}.txt"
                result_file_path = os.path.join(app.config['DOWNLOAD_FOLDER'], result_filename)

                if os.path.exists(result_file_path):
                    result = result_filename
                    print(f"Result file exists at: {result_file_path}")  # Debugging statement

                    subject = 'Your Analysis Results'
                    recipients = [user_email]
                    body = 'Thank you for using our service! Here are the results attached.'

                    msg = MIMEMultipart()
                    msg['Subject'] = subject
                    msg['From'] = app.config['MAIL_DEFAULT_SENDER']
                    msg['To'] = ', '.join(recipients)
                    msg.attach(MIMEText(body, 'plain'))

                    # Attach the result file
                    with open(result_file_path, 'rb') as attachment:
                        part = MIMEBase('application', 'octet-stream')
                        part.set_payload(attachment.read())
                        encoders.encode_base64(part)
                        part.add_header('Content-Disposition', f'attachment; filename="{result_filename}"')
                        msg.attach(part)

                    with smtplib.SMTP(app.config['MAIL_SERVER'], app.config['MAIL_PORT']) as server:
                        server.starttls()
                        server.login(app.config['MAIL_USERNAME'], app.config['MAIL_PASSWORD'])
                        server.sendmail(app.config['MAIL_DEFAULT_SENDER'], recipients, msg.as_string())
                    os.remove(file_path)
                    os.remove(result_file_path)
                else:
                    result = None
                    print("Result file does not exist")
            else:
                flash("Please select a model and provide an email address.", "error")

            if ip in uploaded_files_by_ip:
                uploaded_files_by_ip[ip] += 1
            else:
                uploaded_files_by_ip[ip] = 1

    return render_template('index.html', form=form, result=result)
if __name__ == '__main__':
    app.run(debug=True)
